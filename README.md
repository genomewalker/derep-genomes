
# DeRepG


[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/derep-genomes?include_prereleases&label=version)](https://github.com/genomewalker/derep-genomes/releases) [![derep-genomes](https://github.com/genomewalker/derep-genomes/workflows/derepG_ci/badge.svg)](https://github.com/genomewalker/derep-genomes/actions) [![PyPI](https://img.shields.io/pypi/v/derep-genomes)](https://pypi.org/project/derep-genomes/) [![Conda](https://img.shields.io/conda/v/genomewalker/derep-genomes)](https://anaconda.org/genomewalker/derep-genomes)


DeRepG is a simple genome de-replication tool using [**fastANI**](https://github.com/ParBLiSS/FastANI) and a bit of graph-based analysis. The main objective of DeRepG is to obtain a collection of assemblies that can be used by [**Kraken2**](https://github.com/DerrickWood/kraken2), [**Kaiju**](https://github.com/bioinformatics-centre/kaiju) or [**MMseqs2 taxonomy**](https://github.com/soedinglab/MMseqs2/wiki#taxonomy-assignment) module. It uses a similar approach to the one described in [**Méric, Wick et al. (2019)**](https://www.biorxiv.org/content/10.1101/712166v1) and in [**Parks et al. (2020)**](https://rdcu.be/b3OI7) where they cluster assemblies based on a similarity measure combined with the clustering of the graph after applying a filtering cutoff. In our case, we use as similarity metric the ANI similarity weighted by the fraction aligned, and we don't apply a hard threshold. We follow a different approach where we trim the full connected graph of assemblies per species group, using a re-iterated targeted attack by identifying the edge with the lowest weight that would disconnect the graph if removed. With this approach we let natural patterns emerge in the graph. Afterward, we use the [**Leiden community detection algorithm**](https://www.nature.com/articles/s41598-019-41695-z) to delineate communities in the graph, and then identify the most central assemblies by [**eigenvector centrality**](https://en.wikipedia.org/wiki/Eigenvector_centrality). Then, we extract subgraphs based on the neighbors of those central nodes and we keep those assemblies that are too different (z-score > 2) in terms of genome length and ANI similarity. To speed up the process we do a first MASH distance-based dereplication, where we keep all the edges with a MASH distance <= 0.01. We apply a similar procedure than the one previously described to identify communities and representatives in the MASH graph that will be used in the ANI clustering step. 

While DeRepG can be used to dereplicate a set of genomes and obtain their representative, we cannot assure that it might be the best representative as we are only interested on having an enriched set of genomes for downstream taxonomic classification. For this purpose we recommend to use other approaches like [**Galah**](https://github.com/wwood/galah), [**dRep**](https://drep.readthedocs.io/) or the one described in [**Parks et al. (2020)**](https://rdcu.be/b3OI7).

## Installation

We recommend to have [**conda**](https://docs.conda.io/en/latest/) installed to manage the virtual environments

### Using pip

First we create a conda virtual environment with:

```bash
wget https://raw.githubusercontent.com/genomewalker/derep-genomes/master/environment.yml
conda env create -f environment.yml
```

Then we proceed to install using pip:

```bash
pip install derep-genomes
```

### Using conda

```bash
conda install -c conda-forge -c bioconda -c genomewalker derep-genomes
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
conda env create -f environment.yml
conda activate derep-genomes
pip install -e .
```


## Usage

DeRepG needs a folder containing assemblies in FASTA format and a file describing assemblies with the [**same format**](https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_taxonomy_r95.tsv) as the one distributed by [**GTDB**](/vol/cloud/christiane/NCBI_taxdb_integration/Testing4/viral/assembly_taxonomy.txt). Those are the only requirements. For a full list of option

```
$ derepG --help

usage: derepG [-h] --in-dir DIR --taxa FILE [--db DB] [--threads INT] [--tmp DIR] [--threshold FLOAT] [--slurm-config FILE] [--chunk-size CHUNKS] [--slurm-arr-size INT] [--selected-taxa FILE] [--copy]
              [--out-dir OUT_DIR] [--debug] [--version]

Cluster assemblies in each taxon

required arguments:
  --in-dir DIR          Directory containing all assemblies
  --taxa FILE           TSV file with the taxonomic information

optional arguments:
  -h, --help              show this help message and exit
  --db DB                 SQLite3 DB to store the results (default: derep-genomes.db)
  --threads INT           Number of threads (for fastANI) (default: 16)
  --tmp DIR               Temporary directory (default: /tmp)
  --threshold FLOAT       Z-score filtering threshold (default: 2.0)
  --mash-threhsold FLOAT  Mash distance threshold where to filter (default: 0.01)
  --selected-taxa FILE    File with selected taxa. Taxa should have the following formar: d__;p__;c__;o__;f__;g__;s__
  --copy                  Copy assembly files to the output folder (default: False)
  --out-dir OUT_DIR       Directory where dereplicated assemblies will be copied
  --debug                 Print debug messages (default: False)
  --version               Print program version

SLURM arguments:
  --slurm-config FILE   YAML configuration file for SLURM
  --chunk-size INT      Number of genomes in each chunk for fastANI in SLURM (default: 5)
  --slurm-arr-size INT  Slurm maximum job array size (default: 1000)
  ```

DeRepG uses an SQLite3 database to store the results and the different runs, so if it fails or you have a new set of genomes it doesn't need to rerun the whole dereplication.

One limitation of DeRepG is that fastANI might not be fast when one wants to do all-vs-all comparisons in species with a large number of assemblies. In this case, if you don't want to wait for a long time and you have access to a computer cluster using SLURM it will distribute the fastANI comparisons over the cluster. Check the [**SLURM**](#using-derepg-with-slurm) section for more information

One would run DeRepG as:

```bash
derepG --in-dir data/genomes --taxa data/bac120_taxonomy_r95-1K.tsv --threads 32 --tmp ./ --db test5k-1.db --copy --out-dir gtdb-derep-1k --slurm-config slurm.yaml
```

*--in-dir*: Here we specify the location of the genomes in FASTA format (accepts gzipped files)

*--taxa*: Ataxa file needs to have two columns like the one distributed by GTDB:
> d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica

*--threads*: This means that fastANI runs will be distributed in 16 jobs using 2 threads each

*--tmp*: Location of the temporary folder

*--db*: DeRepG SQLite database location

*--copy* and *--out-dir*: By default DeRepG will generate a TSV file with the kept assemblies and their location. By using *--copy* it will copy the files to the specified folder in *--out-dir*

*--slurm-config*: Location of the SLURM configuration file. Check the [**SLURM**](#using-derepg-with-slurm) section for more info

If everything goes well, you will have your files in the *--out-dir* and a summary file with the columns **taxon**, **accession**, **representative**, **src** and **dst**. If **--copy** has not been set, one can use this file to copy the assembly files a posteriori.

## Using DeRepG with SLURM

Computing the similarities between a large number of genomes can be a time-consuming process. To speed up this step, DeRepG uses [**simple-slurm**](https://github.com/amq92/simple-slurm) to distribute the jobs over a computer cluster if available. By default, DeRepG generates all-vs-all fastANI commands, splits them in groups of n-chunks (**--chunk-size**) and distributes them as array jobs. An example of configuration file can be found [**here**](https://github.com/amq92/simple-slurm#using-configuration-files)
