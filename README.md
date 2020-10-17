
# derep_genomes


[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/derep-genomes?include_prereleases&label=version)](https://github.com/genomewalker/derep-genomes/releases) [![derep-genomes](https://github.com/genomewalker/derep-genomes/workflows/derepG_ci/badge.svg)](https://github.com/genomewalker/derep-genomes/actions) 


A simple genome de-replication tool with fastANI

## Installation

### With [conda](https://docs.conda.io/en/latest/) (recommended)

```bash
conda install -c conda-forge -c bioconda -c genomewalker derep-genomes
```

### With pip

```bash
pip install derep-genomes
```

### Install from source to use the development version

Using pip

```bash
pip install git+ssh://git@github.com/genomewalker/derep-genomes.git@dev
```

By cloning in a dedicated conda environment

```bash
git clone git@github.com:genomewalker/derep-genomes.git
cd derep-genomes
git checkout dev
conda env create -f environment.yml
conda activate derep-genomes
pip install -e .
```
