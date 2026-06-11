import random
from functools import cache
from math import floor, ceil
from pathlib import Path

import numpy as np
import pytest
from conda.testing.helpers import expected_error_prefix

pytest.importorskip("pyteomics")

from workflow.scripts.sample_peptides import (
    cleave_protein_sequence,
    convert_fasta_str_to_dict,
    sample_accessions,
    sample_noise_peptides,
    sample_peptides,
    sample_peptides_from_sequence,
    UniProtConnector,
)


class MockedUniprotConnector:
    def __init__(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def _generate_random_peptide_sequence(self) -> str:
        amino_acid_alphabet = "ARNDCEQGHILKMFPSTWYV"
        length = random.randint(10, 2*10**4)
        amino_acids = random.choices(amino_acid_alphabet, k=length)
        return "".join(amino_acids)

    @cache
    def _get_sequence_for_acc(self, acc) -> str:
        return self._generate_random_peptide_sequence()

    def get_accession_to_sequence_mapping(self, noise_accessions: list[str]):
        acc2seq: dict[str, str] = {}
        for acc in noise_accessions:
            acc2seq[acc] = self._get_sequence_for_acc(acc)
        return acc2seq


def test_get_fasta() -> None:
    with UniProtConnector() as connector:
        fasta = connector.get_fasta(["A0A3R6E0N5"])
    assert fasta == """\
>tr|A0A3R6E0N5|A0A3R6E0N5_9FIRM Glycosyl transferase OS=Roseburia intestinalis OX=166486 GN=DW264_18045 PE=4 SV=1
MCGGDDILKYRTYCKNQRDVAFVINGIIDEYWCGKLSEKEMKEDILTLYENNKEKLFKDG
QFTKIIQQQCGKKRINVISQILKNKLEKLE
"""


def test_get_fasta_returns_more_than_25_results() -> None:
    accessions = [
        'A0A3R6E0N5', 'A0A3E4LLW1', 'A0A7L6WN97', 'A0A6S6M774', 'A0A174UR68', 'A0A3R9NJY9', 'A0A0B0UEA8', 'A0A5C1JJV1',
        'A0A3E5GV22', 'A0A0L0LSS8', 'A0A427ZU75', 'A0A6H9PP43', 'A0A6G6L3K2', 'A0A6N3APR0', 'A0A0A1GRF5', 'A0A413EJW6',
        'A0A7I9AKF9', 'A0A1Q6B4J7', 'A0A558LVQ5', 'A0A2N5PYJ1', 'A0A139L2U9', 'F7LRI5', 'A0A376TL63', 'A0A8T3LEX4',
        'A0A139L7B9', 'A0A1S6GKC3', 'A0A7M1NVS6', 'A0A7D7DSC8', 'A0A0M1UJ95', 'A0A174PW89', 'A0A448PJ07', 'A0A415D293',
        'A0A7J4XQ64', 'A0A412P7C0', 'A0A0G3B8Y9', 'A0A9Q6F4E3', 'A0A7W9SEH0', 'A0A137SRT4', 'A0A2H4U0C6', 'A0A412RPZ4',
    ]
    with UniProtConnector() as connector:
        fasta = connector.get_fasta(accessions)
    assert fasta.count(">") > 25


def test_sample_peptides_from_sequence():
    expected = ['TGPNLHGLFGR', 'CSQCHTVEK', 'GIIWGEDTLME', 'TGQAPGYSYTAANK', 'IFIMKCSQCHTVEK']
    random.seed(12345)
    numpy_random_number_generator = np.random.default_rng(seed=random.randint(0, 2**31-1))
    sequence = "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLME"
    returned = sample_peptides_from_sequence(numpy_random_number_generator, sequence)
    assert returned == expected


@pytest.mark.parametrize(
    ("missed_cleavages", "expected"),
    [
        (0, {'TGQAPGYSYTAANK', 'TGPNLHGLFGR', 'GIIWGEDTLME', 'CSQCHTVEK'}),
        (1, {'CSQCHTVEKGGK', 'KTGQAPGYSYTAANK', 'HKTGPNLHGLFGR', 'TGQAPGYSYTAANKNK', 'TGPNLHGLFGRK', 'IFIMKCSQCHTVEK', 'NKGIIWGEDTLME', 'MGDVEKGK'}),
        (2, {'KTGQAPGYSYTAANKNK', 'GKKIFIMK', 'TGPNLHGLFGRKTGQAPGYSYTAANK', 'TGQAPGYSYTAANKNKGIIWGEDTLME', 'KIFIMKCSQCHTVEK', 'IFIMKCSQCHTVEKGGK', 'GGKHKTGPNLHGLFGR', 'HKTGPNLHGLFGRK', 'MGDVEKGKK', 'CSQCHTVEKGGKHK'}),
    ]
)
def test_cleave_protein_sequence(missed_cleavages: int, expected: set[str]):
    sequence = "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLME"
    assert cleave_protein_sequence(sequence, missed_cleavages=missed_cleavages) == expected


def test_convert_fasta_str_to_dict():
    fasta = """>sp|P12345|AATM_RABIT Aspartate aminotransferase, mitochondrial OS=Oryctolagus cuniculus OX=9986 GN=GOT2 PE=1 SV=2
MALLHSARVLSGVASAFHPGLAAAASARASSWWAHVEMGPPDPILGVTEAYK
>sp|P99999|CYC_HUMAN Cytochrome c OS=Homo sapiens OX=9606 GN=CYCS PE=1 SV=2
MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIW
GEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE"""
    assert convert_fasta_str_to_dict(fasta) == {
        "P12345": "MALLHSARVLSGVASAFHPGLAAAASARASSWWAHVEMGPPDPILGVTEAYK",
        "P99999": "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE",
    }


def test_sample_peptides(workflow_path: Path, tmp_path: Path) -> None:
    expected_noise_percentage = 10
    test_data = Path(__file__).parent.parent / "rules" / "simulate_sample/sample_peptides/data/results/simulation/"

    no_noise = tmp_path / "no_noise.txt"
    sample_peptides(test_data / "tax2accessions.json", test_data / "sample_taxons_lineage_1.tsv", no_noise, noise_percentage=0)
    signal_peptides = [line.strip() for line in no_noise.read_text(encoding="utf-8").splitlines()]

    with_noise = tmp_path / "with_noise.txt"
    sample_peptides(test_data / "tax2accessions.json", test_data / "sample_taxons_lineage_1.tsv", with_noise, noise_percentage=expected_noise_percentage)
    signal_and_noise_peptides = [line.strip() for line in with_noise.read_text(encoding="utf-8").splitlines()]

    signal_set = set(signal_peptides)
    signal_and_noise_set = set(signal_and_noise_peptides)
    returned_noise_percentage = ceil((len(signal_and_noise_peptides) - len(signal_peptides)) / len(signal_peptides) * 100)

    assert signal_set.issubset(signal_and_noise_set)
    assert len(signal_and_noise_set - signal_set) >= 1
    assert returned_noise_percentage == expected_noise_percentage


def test_sample_noise_peptides():
    random.seed(123)
    numpy_random_number_generator = np.random.default_rng(seed=random.randint(0, 2 ** 31 - 1))

    tax2acc = {"123": ["foo"]}
    noise_taxid_population = list(tax2acc)
    signal_peptides = ["IGHDNFRK"]

    with MockedUniprotConnector() as connector:
        peptides = sample_noise_peptides(
            26,
            noise_taxid_population,
            tax2acc,
            signal_peptides,
            numpy_random_number_generator,
            connector
        )
    assert len(peptides) == 26
    assert signal_peptides[0] not in peptides


def test_sample_noise_peptides_raises_runtime_error(monkeypatch: pytest.MonkeyPatch) -> None:
    random.seed(123)
    numpy_random_number_generator = np.random.default_rng(seed=random.randint(0, 2 ** 31 - 1))

    tax2acc = {"123": ["foo"]}
    noise_taxid_population = list(tax2acc)
    sequence_without_cleavage_site: str = "ACDEFGHILMNPQSTVWYV"
    signal_peptides = []

    # Define replacement function
    def fake_method(self: MockedUniprotConnector, accessions: list[str]) -> dict[str, str]:
        return {acc: sequence_without_cleavage_site for acc in accessions}

    # Apply monkeypatch
    monkeypatch.setattr(
        MockedUniprotConnector,
        "get_accession_to_sequence_mapping",
        fake_method,
    )

    with MockedUniprotConnector() as connector:
        with pytest.raises(RuntimeError, match="Maximum number of tries reached: 10"):
            sample_noise_peptides(
                2,
                noise_taxid_population,
                tax2acc,
                signal_peptides,
                numpy_random_number_generator,
                connector
            )



def test_sample_accessions() -> None:
    tax_ids = ["1", "2", "3"]
    tax2acc = {
        "1": ["1a", "1b", "1c"],
        "2": ["2a", "2b", "2c"],
        "3": ["3a", "3b", "3c"],
    }
    random.seed(12345)
    returned = sample_accessions(tax2acc, tax_ids)
    assert returned == ['1b', '1a', '1c', '2b', '2a', '2c', '3c', '3b', '3a']


def test_get_accession_to_sequence_mapping():
    accessions = [
        'A0A3R6E0N5', 'A0A3E4LLW1', 'A0A7L6WN97', 'A0A6S6M774', 'A0A174UR68', 'A0A3R9NJY9', 'A0A0B0UEA8', 'A0A5C1JJV1',
        'A0A3E5GV22', 'A0A0L0LSS8', 'A0A427ZU75', 'A0A6H9PP43', 'A0A6G6L3K2', 'A0A6N3APR0', 'A0A0A1GRF5', 'A0A413EJW6',
        'A0A7I9AKF9', 'A0A1Q6B4J7', 'A0A558LVQ5', 'A0A2N5PYJ1', 'A0A139L2U9', 'F7LRI5', 'A0A376TL63', 'A0A8T3LEX4',
        'A0A139L7B9', 'A0A1S6GKC3', 'A0A7M1NVS6', 'A0A7D7DSC8', 'A0A0M1UJ95', 'A0A174PW89', 'A0A448PJ07', 'A0A415D293',
        'A0A7J4XQ64', 'A0A412P7C0', 'A0A0G3B8Y9', 'A0A9Q6F4E3', 'A0A7W9SEH0', 'A0A137SRT4', 'A0A2H4U0C6', 'A0A412RPZ4',
        'A0A8G1S9B2', 'A0A413I9L1', 'A0A2X6D6G0', 'O32560', 'A0A139KIL5', 'A0A9Q8DRL9', 'A0A250KGK4', 'A0A377DKM1',
        'A0A2N0UJJ5', 'A0A772E7H1', 'A0A3E5EFK6', 'A0A930LQ42', 'A0A827DYC8', 'A0A9P2HL15', 'A0A2N0URT8', 'A0A9Q4IVR9',
        'A0A4R6CTN2', 'A0A6P1Y4N3', 'A0A8B3B8D5', 'A0A3R5ZQX8', 'A0A174ERA0', 'A0A6D0FWW0', 'J1LM39', 'A0A7M1NXV3',
        'A0A6A1Z4W4', 'A0A3S4K6R1', 'A0A1V3CDT4', 'A0A4Q5HEZ6', 'A0A8T5ZQ29', 'A0A3E5ET67', 'A0A4Q5IDM6', 'A0A413YX76',
        'A0A1V2FXF9', 'A0A5M6A476', 'A0A2U2RUX8', 'A0A849YN47', 'A0A2A3UQ42', 'A0A0L0LST0', 'A0A828ABF2', 'A0A826J2R3',
        'A0A376K0H5', 'A0A174NBZ0', 'A0A139LVH8', 'A0A930Q174', 'A0A2N0SPF6', 'A0A1Q8I2P6', 'A0A2X1JGP1', 'A0A2X3K270',
        'A0A7J4XUU1', 'A0A173RK00', 'A0A827QNQ9', 'A0A9Q7ZIA1', 'A0A378TZE8', 'A0A921K5C8'
    ]
    with UniProtConnector() as uniprot_connector:
        accession_to_sequence_mapping = uniprot_connector.get_accession_to_sequence_mapping(accessions)
    assert len(accession_to_sequence_mapping) > 25
    # assert set(accessions) == set(accession_to_sequence_mapping)
