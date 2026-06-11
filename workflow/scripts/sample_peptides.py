import random
from math import floor
from typing import TYPE_CHECKING, Generator, Self
from hashlib import sha256

import numpy as np
import pandas as pd
from itertools import groupby
from typing import cast
from pathlib import Path
from numpy.random import Generator as NumpyGenerator
from pyteomics.parser import cleave
import requests

import json
import sys

if TYPE_CHECKING:
    from snakemake import script

LOG_HANDLE = sys.stderr


class UniProtConnector:
    url = "https://rest.uniprot.org/uniprotkb/"

    def __init__(self):
        self.session = requests.Session()
        self.uniprot_version = requests.get(self.url).headers.get("X-UniProt-Release")
        print(f"Uniprot Release Number: {self.uniprot_version}", file=LOG_HANDLE)
        self._accession_to_sequence_cache: dict[str, str] = {}

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.session.close()

    def get_fasta(self, uniprot_accessions: list[str]) -> str:
        fastas = []
        chunk_size = min(25, len(uniprot_accessions))
        for chunk_start in range(0, len(uniprot_accessions), chunk_size):
            chunk = uniprot_accessions[chunk_start: chunk_start + chunk_size]
            accession_string = ["accession%3A" + acc for acc in chunk]
            url = self.url + f"search?query={'+OR+'.join(accession_string)}&format=fasta"
            response = self.session.get(url)
            response.raise_for_status()
            fastas.append(response.text)
        return "".join(fastas)

    def get_accession_to_sequence_mapping(self, accessions: list[str]) -> dict[str, str]:
        missing_accessions = list(set(accessions) - set(self._accession_to_sequence_cache))
        fasta = self.get_fasta(missing_accessions)
        retrieved_accession_to_sequence_mapping = convert_fasta_str_to_dict(fasta)
        unretrievable_accessions = sorted({acc for acc in missing_accessions if acc not in retrieved_accession_to_sequence_mapping})
        print(f"Could not retrieve sequences for the following accessions: {' '.join(unretrievable_accessions)}", file=LOG_HANDLE)
        self._accession_to_sequence_cache.update(retrieved_accession_to_sequence_mapping)
        return {
            acc: self._accession_to_sequence_cache[acc]
            for acc in accessions
            if acc not in unretrievable_accessions
        }


def sample_accessions(tax2acc: dict[str, list[str]], tax_ids: list[str]) -> list[str]:
    """Sample accessions for each tax_id.

    Samples min(len(accessions_of_this_tax_id), 100) accessions for each tax id in tax_ids.

    Args:
        tax2acc: mapping of taxonomic identifiers to list of accessions
        tax_ids: list of taxonomic identifiers to sample accessions

    Returns:
        list of sampled accessions
    """
    accessions_sample: list[str] = []
    for tax_id in tax_ids:
        accessions = tax2acc[tax_id]
        sample_size_for_current_tax_id = min(len(accessions), 100)
        accessions_sample.extend(random.sample(accessions, sample_size_for_current_tax_id))
    return accessions_sample


def sample_peptides(
    tax2acc_map_file: str | Path,
    lineage_file: str | Path,
    output: str | Path,
    noise_percentage: int = 1
):
    """Generate sample of peptides for tax_ids in the provided lineage file.

    The random module is seeded with the lineage_file.
    Adds false positive peptides to the end of the output.

    Args:
        tax2acc_map_file: path to json with tax_ids as keys and list of accessions as values
        lineage_file: path to tab seperated lineage file with column "query", which is used to look up accessions to
          sample from in tax2acc_map_file and as randomness seed
        output: path to output file, containing one peptide per line
        noise_percentage: percentage of the number of sampled signal peptides to add to the sample
    """
    with open(tax2acc_map_file, "r") as handle:
        tax2acc = cast(dict[str, list[str]], json.load(handle))
    df_tax_ids = pd.read_csv(lineage_file, usecols=["query"], dtype={"query": str}, sep="\t")
    signal_tax_ids: list[str] = df_tax_ids["query"].to_list()

    # sample accessions
    salted_seed = str(lineage_file) + "sample_peptides"
    seed = sha256(salted_seed.encode("utf-8")).hexdigest()
    random.seed(seed)
    numpy_random_number_generator = np.random.default_rng(seed=random.randint(0, 2 ** 31 - 1))
    with UniProtConnector() as uniprot_connector:
        peptides = sample_signal_peptides(signal_tax_ids, tax2acc, numpy_random_number_generator, uniprot_connector)

        if noise_percentage:
            noise_tax_id_population = [
                tax_id
                for tax_id in tax2acc
                if tax_id not in signal_tax_ids
            ]
            noise_peptide_count = floor(noise_percentage / 100 * len(peptides))
            peptides.extend(
                sample_noise_peptides(
                    noise_peptide_count=noise_peptide_count,
                    noise_tax_id_population=noise_tax_id_population,
                    tax2acc=tax2acc,
                    signal_peptides=peptides,
                    numpy_random_number_generator=numpy_random_number_generator,
                    uniprot_connector=uniprot_connector
                )
            )

    # write output
    with open(output, "w") as handle:
        handle.write("\n".join(peptides) + "\n")


def sample_signal_peptides(
    signal_tax_id_population: list[str],
    tax2acc: dict[str, list[str]],
    numpy_random_number_generator: NumpyGenerator,
    uniprot_connector: UniProtConnector
) -> list[str]:
    accessions_sample = sample_accessions(tax2acc, signal_tax_id_population)

    print(f"Sampled {len(accessions_sample)} signal accessions", file=LOG_HANDLE)

    acc2seq = uniprot_connector.get_accession_to_sequence_mapping(accessions_sample)

    # cleave proteins sequences and sample peptides
    sequences = sorted(acc2seq.values())
    peptides = list()
    for seq in sequences:
        peptides.extend(
            sample_peptides_from_sequence(numpy_random_number_generator, seq)
        )
    print(f"Sampled {len(peptides)} signal peptides", file=LOG_HANDLE)
    return peptides


def sample_noise_peptides(
    noise_peptide_count: int,
    noise_tax_id_population: list[str],
    tax2acc: dict[str, list[str]],
    signal_peptides:  list[str],
    numpy_random_number_generator: NumpyGenerator,
    uniprot_connector: UniProtConnector
) -> list[str]:
    accession_population = [
        acc
        for tax_id in noise_tax_id_population
        for acc in tax2acc[tax_id]
    ]

    def random_sequence_generator() -> Generator[str, None, None]:
        while True:
            noise_accessions = random.choices(accession_population, k=25)
            noise_acc2seq = uniprot_connector.get_accession_to_sequence_mapping(noise_accessions)
            yield from noise_acc2seq.values()

    # cleave proteins sequences and sample peptides
    maximum_number_of_tries = 5*noise_peptide_count
    noise_peptides = list()
    peptides_drawn = 0
    for executed_number_of_tries, seq in enumerate(random_sequence_generator()):
        if executed_number_of_tries >= maximum_number_of_tries:
            raise RuntimeError(f"Maximum number of tries reached: {maximum_number_of_tries}")
        peptide_sample_current_sequence = sample_peptides_from_sequence(numpy_random_number_generator, seq)
        if not peptide_sample_current_sequence:
            continue
        peptide = random.choice(
                peptide_sample_current_sequence
            )
        if peptide in noise_peptides:
            continue
        if peptide in signal_peptides:
            continue
        noise_peptides.append(peptide)
        peptides_drawn += 1
        if peptides_drawn == noise_peptide_count:
            break
    print(f"Sampled {len(noise_peptides)} noise peptides", file=LOG_HANDLE)
    return noise_peptides


def test_sample_peptides(tmpdir: Path) -> None:
    test_data_path = Path(
        __file__).parent.parent.parent / "test/unit/simulate_sample/sample_peptides/data/results/simulation"
    test_data_path = Path("/home/hennings/Downloads/delete_me/")
    sample_peptides(
        tax2acc_map_file=test_data_path / "tax2accessions.json",
        lineage_file=test_data_path / "sample_taxons_lineage_0.tsv",
        output=tmpdir / "peptides.txt",
    )


def sample_peptides_from_sequence(
    numpy_random_number_generator: NumpyGenerator, sequence: str
) -> list[str]:
    """Digest sequence according to tryptic digestion rules and return a random sample of unique peptides.

    The sample size is drawn from a poisson distribution with mean 10.
    Digestion allows up to 2 missed cleavages, with the weights of each missed cleavage as follows:
    - 0 missed cleavages: 100
    - 1: 10
    - 2: 1

    Args:
        numpy_random_number_generator: initialized numpy Generator
        sequence: amino acid sequence to be digested

    Returns:
        sample of peptides
    """
    (
        peptides_0_missed_cleavages,
        peptides_1_missed_cleavages,
        peptides_2_missed_cleavages,
    ) = [
        list(
            sorted(cleave_protein_sequence(sequence, missed_cleavages=missed_cleavages))
        )
        for missed_cleavages in [0, 1, 2]
    ]
    sample_size = numpy_random_number_generator.poisson(lam=10)
    missed_cleavages_distribution = random.choices(
        [0, 1, 2], weights=[100, 10, 1], k=sample_size
    )

    return (
        random.sample(
            population=peptides_0_missed_cleavages,
            k=min(
                missed_cleavages_distribution.count(0), len(peptides_0_missed_cleavages)
            ),
        )
        + random.sample(
            population=peptides_1_missed_cleavages,
            k=min(
                missed_cleavages_distribution.count(1), len(peptides_1_missed_cleavages)
            ),
        )
        + random.sample(
            population=peptides_2_missed_cleavages,
            k=min(
                missed_cleavages_distribution.count(2), len(peptides_2_missed_cleavages)
            ),
        )
    )


def cleave_protein_sequence(
    sequence: str,
    min_length: int = 7,
    max_length: int = 30,
    missed_cleavages: int = 0,
    *args,
    **kwargs,
) -> set[str]:
    kwargs.update(
        {
            "sequence": sequence,
            "rule": "trypsin",
            "min_length": min_length,
            "max_length": max_length,
        }
    )
    if missed_cleavages > 0:
        return cleave(*args, missed_cleavages=missed_cleavages, **kwargs) - cleave(
            *args, missed_cleavages=missed_cleavages - 1, **kwargs
        )
    return cleave(*args, missed_cleavages=missed_cleavages, **kwargs)


def convert_fasta_str_to_dict(fasta):
    protein_dict = {}
    fasta = fasta.strip()
    fasta_iter = (x[1] for x in groupby(fasta.split("\n"), lambda line: line[0] == ">"))

    for header in fasta_iter:
        # keep only uniprot accession, between pipe-chars ('|')
        acc = header.__next__().split("|")[1]

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in fasta_iter.__next__())

        protein_dict[acc] = seq
    return protein_dict


snakemake: "script.Snakemake"
if snakemake := globals().get("snakemake"):
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        sample_peptides(
            tax2acc_map_file=snakemake.input["tax2acc_map"],
            lineage_file=snakemake.input["lineage"],
            output=snakemake.output[0],
            noise_percentage=int(snakemake.wildcards["percentage"]),
        )
