import pandas as pd
import pyfastx
import os
from derep_genomes.general import (
    find_all_assemblies,
    load_classifications,
    shorten_accession,
    get_accession,
)
from derep_genomes import __version__
import logging

log = logging.getLogger("my_logger")
log.setLevel(logging.INFO)


def rename_header(header):
    pass


def read_fasta(data):
    print("hello")
    file = str(data[2])
    print(file)
    fa = pyfastx.Fasta()
    print(fa)


def index(args):
    all_assemblies = find_all_assemblies(args.in_dir)
    classifications = load_classifications(args.tax_file)

    log.info("Filtering taxa with assemblies")
    classifications_df = pd.DataFrame(
        [
            (key, var, shorten_accession(var))
            for (key, L) in classifications.items()
            for var in L
        ],
        columns=["taxon", "accession", "accession_nover"],
    )
    log.info(
        "Found {:,} assembly files in {}".format(
            classifications_df.shape[0], args.in_dir
        )
    )
    all_assemblies_df = list(map(get_accession, all_assemblies))
    assm_data = classifications_df.merge(pd.DataFrame(all_assemblies_df))
    assm_data.loc[:, "file"] = assm_data["assembly"].apply(os.path.basename)

    assm_files = list(
        assm_data[["accession", "accession_nover", "assembly"]].to_records(index=False)
    )
    read_fasta(assm_files[0])
    map(read_fasta, assm_files)
