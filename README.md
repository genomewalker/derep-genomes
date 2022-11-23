
# DeRepG


[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/derep-genomes?include_prereleases&label=version)](https://github.com/genomewalker/derep-genomes/releases) [![derep-genomes](https://github.com/genomewalker/derep-genomes/workflows/derepG_ci/badge.svg)](https://github.com/genomewalker/derep-genomes/actions) [![PyPI](https://img.shields.io/pypi/v/derep-genomes)](https://pypi.org/project/derep-genomes/)

DeRepG is a simple genome de-replication tool using [**fastANI**](https://github.com/ParBLiSS/FastANI) and a bit of graph-based analysis. The main objective of DeRepG is to obtain a collection of assemblies that can be used by [**Kraken2**](https://github.com/DerrickWood/kraken2), [**Kaiju**](https://github.com/bioinformatics-centre/kaiju) or [**MMseqs2 taxonomy**](https://github.com/soedinglab/MMseqs2/wiki#taxonomy-assignment) module. It uses a similar approach to the one described in [**Méric, Wick et al. (2019)**](https://www.biorxiv.org/content/10.1101/712166v1) and in [**Parks et al. (2020)**](https://rdcu.be/b3OI7) where they cluster assemblies based on a similarity measure combined with the clustering of the graph after applying a filtering cutoff. In our case, we use the ANI similarity weighted by the fraction aligned as a metric, and we don't apply a hard threshold. We follow a different approach where we trim the full connected graph of assemblies per species group, using a re-iterated targeted attack by identifying the edge with the lowest weight that would disconnect the graph if removed. With this approach we let natural patterns emerge in the graph. Afterward, we use the [**Leiden community detection algorithm**](https://www.nature.com/articles/s41598-019-41695-z) to delineate communities in the graph, and then identify the most central assemblies by [**eigenvector centrality**](https://en.wikipedia.org/wiki/Eigenvector_centrality). Then, we extract subgraphs based on the neighbors of those central nodes and we keep those assemblies that are too different (z-score > 2) in terms of genome length and ANI similarity. To speed up the process we do a first MASH or Dashing distance-based dereplication, where we keep all the edges with a distance <= 0.01. We apply a similar procedure to the one previously described to identify communities and representatives in the MASH/Dashing graph that will be used in the ANI clustering step. 

While DeRepG can be used to dereplicate a set of genomes and obtain their representative, we cannot assure that it might be the best representative as we are only interested in having an enriched set of genomes for downstream taxonomic classification. For this purpose, we recommend using other approaches like [**Galah**](https://github.com/wwood/galah), [**dRep**](https://drep.readthedocs.io/) or the one described in [**Parks et al. (2020)**](https://rdcu.be/b3OI7).

> As fastANI produces [asymmetric similarity scores](https://github.com/ParBLiSS/FastANI/issues/36) when doing all-vs-all comparisons, we use the mean of the two scores as the similarity score and we select the largest aligned fraction for the graph analysis.

## Installation
We recommend having **conda** installed to manage the virtual environments

### Using pip

First, we create a conda virtual environment with:

```bash
wget https://raw.githubusercontent.com/genomewalker/derep-genomes/master/environment.yml
conda env create -f environment.yml
```

Then we proceed to install using pip:

```bash
pip install derep-genomes
```

### Install from source to use the development version

Using pip

```bash
pip install git+ssh://git@github.com/genomewalker/derep-genomes.git
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

usage: derepG [-h] --data FILE [--db DB] [--workers INT] [--threads INT] [--prefix PREFIX] [--tmp DIR] [--threshold FLOAT] [--xash-threshold FLOAT] [--max-assemblies INT] [--slurm-config FILE] [--chunk-size INT]
              [--slurm-arr-size INT] [--selected-taxa FILE] [--copy] [--dashing] [--out-dir OUT_DIR] [--min-genome-size INT] [--ani-fraglen-fraction FLOAT] [--debug] [--version]

Cluster assemblies in each taxon

required arguments:
  --data FILE           TSV file with the genome information

optional arguments:
  -h, --help            show this help message and exit
  --db DB               SQLite3 DB to store the results (default: derep-genomes.db)
  --workers INT         Number of workers to use (default: 2)
  --threads INT         Number of threads to use (default: 2)
  --prefix PREFIX       Prefix for the file name results. If not assigned it uses Ymd-HMS (default: 20221012-124218)
  --tmp DIR             Temporary directory (default: /tmp)
  --threshold FLOAT     Z-score filtering threshold (default: 2.0)
  --xash-threshold FLOAT
                        Mash/Dashing distance threshold where to filter (default: 0.01)
  --max-assemblies INT  Maximum number of assemblies to process for small batches (default: 10)
  --selected-taxa FILE  File with selected taxa. Taxa should have the following formar: d__;p__;c__;o__;f__;g__;s__
  --copy                Copy assembly files to the output folder (default: False)
  --dashing             Use Dashing instead of Mash (default: False)
  --out-dir OUT_DIR     Directory where dereplicated assemblies will be copied
  --min-genome-size INT
                        Minimum genome size where to apply heuristics to find ANI fragment size (default: 20000000)
  --ani-fraglen-fraction FLOAT
                        Fraction of the genome size used to estimate the ANI used fragment length (default: 0.005)
  --debug               Print debug messages (default: False)
  --version             Print program version

SLURM arguments:
  --slurm-config FILE   YAML configuration file for SLURM
  --chunk-size INT      Number of genomes in each chunk for fastANI i
  ```

DeRepG uses an SQLite3 database to store the results and the different runs, so if it fails or you have a new set of genomes it doesn't need to rerun the whole dereplication.

One limitation of DeRepG is that fastANI might not be fast when one wants to do all-vs-all comparisons in species with a large number of assemblies. In this case, if you don't want to wait for a long time and you have access to a computer cluster using SLURM it will distribute the fastANI comparisons over the cluster. Check the [**SLURM**](#using-derepg-with-slurm) section for more information.

One would run DeRepG as:

```bash
derepG --data data/genomes --taxa data/bac120_taxonomy_r95-1K.tsv --workers 2 --threads 32 --tmp ./ --db test5k-1.db --copy --out-dir gtdb-derep-1k --slurm-config slurm.yaml
```

Where:

*--data*: A TSV file where we specify the location of the genomes in FASTA format (accepts gzipped files). The TSV must have the following columns: `accession,taxonomy,file`


```bash
accession       taxonomy        file
GCF_000566285.1 d__Bacteria;l__Bacteria;k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli       /data/genomes/GCF_000566285.1_genomic.fna.gz
GCF_003460375.1 d__Bacteria;l__Bacteria;k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli       /data/genomes/GCF_003460375.1_genomic.fna.gz
```


*--workers*: Number of workers to use. This is the number of taxa to process in parallel. If you have a large number of taxa, you might want to increase this number.

*--threads*: Total number of threads available for the workers. This is the number of threads used by fastANI and Mash/Dashing.

*--tmp*: Location of the temporary folder

*--db*: DeRepG SQLite database location

*--copy* and *--out-dir*: By default DeRepG will generate a TSV file with the kept assemblies and their location. By using *--copy* it will copy the files to the specified folder in *--out-dir*

*--slurm-config*: Location of the SLURM configuration file. Check the [**SLURM**](#using-derepg-with-slurm) section for more info

If everything goes well, you will have your files in the *--out-dir* and a summary file with the columns **taxon**, **accession**, **representative**, **src** and **dst**. If **--copy** has not been set, one can use this file to copy the assembly files a posteriori.

## Using DeRepG with SLURM

Computing the similarities between a large number of genomes can be a time-consuming process. To speed up this step, DeRepG uses [**simple-slurm**](https://github.com/amq92/simple-slurm) to distribute the jobs over a computer cluster if available. By default, DeRepG generates all-vs-all fastANI commands, splits them in groups of n-chunks (**--chunk-size**) and distributes them as array jobs. An example of configuration file can be found [**here**](https://github.com/amq92/simple-slurm#using-configuration-files)
