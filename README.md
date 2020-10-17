
# derep_genomes


[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/derep-genomes?include_prereleases&label=version)](https://github.com/genomewalker/derep-genomes/releases) [![derep-genomes](https://github.com/genomewalker/derep-genomes/workflows/derepG_ci/badge.svg)](https://github.com/genomewalker/derep-genomes/actions) 


A simple genome de-replication tool with fastANI

## Installation

Clone the repo:

```bash
git clone https://github.com/genomewalker/derep-genomes.git
cd derep-genomes
```


Create a conda env:

```bash
conda env create -f environment.yaml
conda activate genome-derep
```

Install it:
```bash
python setup.py install 
```


# PyDamage

Pydamage, is a Python software to automate the process of contig damage identification and estimation.
After modelling the ancient DNA damage using the C to T transitions, Pydamage uses a likelihood ratio test to discriminate between truly ancient, and modern contigs originating from sample contamination.

## Installation

### With [conda](https://docs.conda.io/en/latest/) (recommended)

```bash
conda install -c conda-forge -c bioconda -c maxibor pydamage
```

### With pip

```bash
pip install pydamage
```

### Install from source to use the development version

Using pip

```bash
pip install git+ssh://git@github.com/maxibor/pydamage.git@dev
```

By cloning in a dedicated conda environment

```bash
git clone git@github.com:maxibor/pydamage.git
cd pydamage
git checkout dev
conda env create -f environment.yml
conda activate pydamage
pip install -e .
```


## Quick start

```bash
pydamage aligned.bam reference.fa
```

## CLI help

Command line interface help message

```bash
pydamage --help
```

## Documentation

[pydamage.readthedocs.io](https://pydamage.readthedocs.io)
