# megadudes-evaluation workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.25.0-brightgreen.svg)](https://snakemake.bitbucket.io)

Evaluation workflow for taxonomic profiling of metaproteomic data using the [DUDes](https://github.com/pirovc/dudes) and [Unipept](https://unipept.ugent.be/). The workflow can simulate peptide data and analyze MSFragger output.

## Usage

### Setup

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository.

1. install miniconda according to their [instructions](https://www.anaconda.com/docs/getting-started/miniconda/install/overview)
2. clone this repo
3. install Snakemake: `conda env create -f env.yaml`
4. download databases ⚠️requires ~100gb of disk space⚠️: `./download_resources.sh`

### Configuration

Configure the workflow by editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.tsv` to specify your sample setup.

#### `config.yaml`

| Parameter | Description |
| --- | --- |
| `query_dbs` | List of paths to protein sequence database files (e.g. UniProt `.fasta.gz` files) used as reference databases for sequence alignment and taxonomic profiling with DUDes. |
| `alignment_methods` | List of alignment methods to run against each `query_dbs` entry. Supported values: `diamond`, `mmseqs2`, `mmseqs2_top_10`. |

#### `samples.tsv`

Tab-separated file with one row per sample to be evaluated. Required columns:

| Column | Description |
| --- | --- |
| `sample_name` | Unique identifier for the sample. |
| `ground_truth` | Path to a tsv file describing the expected taxonomic composition of the sample. Must contain the columns `superkingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `subspecies`, each holding the corresponding NCBI taxonomy ID (additional columns are ignored). |
| `msfragger_peptides_tsv` | Path to the `peptide.tsv` output of an MSFragger search. Only the `Peptide` column is used. |

### Execution

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N --resources mem_gb=21

## Updating

- update the local environment: `conda update -n snakemake --all`
- update env file: `conda env export -n snakemake > env.yaml`
- update lock-file: `conda list --explicit -n snakemake > spec-file.txt`
- update pre-commit hooks: `pre-commit autoupdate`

## Contributing

- setup pre-commit hooks: `pre-commit install`
- execute tests: `pytest`
