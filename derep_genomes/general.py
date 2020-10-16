import argparse
import sys
import shutil
import collections
import gzip
import os
import pathlib
import collections
import textwrap
import logging
import pandas as pd
from multiprocessing import Pool
from functools import partial
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from os import devnull
import tqdm
from derep_genomes import __version__


log = logging.getLogger("my_logger")
log.setLevel(logging.INFO)


def is_debug():
    return logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG


# From: https://stackoverflow.com/a/11541450
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        if os.path.isfile(arg):
            parser.error("The file %s does not exist!" % arg)
        else:
            parser.error("The directory %s does not exist!" % arg)
    else:
        if os.path.isfile(arg):
            return open(arg, "r")  # return an open file handle
        else:
            return arg


def get_arguments():

    parser = argparse.ArgumentParser(description="Cluster assemblies in each taxon")
    parser.add_argument(
        "--in-dir",
        dest="in_dir",
        required=True,
        metavar="DIR",
        type=lambda x: is_valid_file(parser, x),
        help="Directory containing all GTDB assemblies",
    )
    parser.add_argument(
        "--tmp",
        type=str,
        default="/tmp",
        dest="tmp_dir",
        help="Tmp directory where dereplicated assemblies will be copied",
    )
    parser.add_argument(
        "--taxa",
        required=True,
        metavar="FILE",
        type=lambda x: is_valid_file(parser, x),
        dest="tax_file",
        help="GTDB taxonomy file",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=2,
        help="Z-score filtering threshold",
    )
    parser.add_argument(
        "-p", "--threads", type=int, default=16, help="Number of threads (for fastANI)"
    )
    parser.add_argument(
        "-c",
        "--chunk-size",
        dest="chunks",
        type=int,
        default=5,
        help="Number of genomes in each chunk for fastANI in SLURM",
    )
    parser.add_argument(
        "-a",
        "--slurm-max-array-size",
        dest="max_jobs_array",
        type=int,
        default=1000,
        help="Slurm maximum job array size",
    )

    parser.add_argument(
        "--db",
        dest="db",
        type=str,
        default="derep-genomes.db",
        help="SQLite3 DB to store the results",
    )
    parser.add_argument(
        "-s",
        "--slurm_config",
        dest="slurm_config",
        required=False,
        help="YAML configuration for slurm",
        metavar="FILE",
        type=lambda x: is_valid_file(parser, x),
    )

    parser.add_argument(
        "--selected-taxa",
        dest="selected_taxa",
        required=False,
        help="File with selected taxa. Taxa should have the following formar: d__;p__;c__;o__;f__;g__;s__",
        metavar="FILE",
        type=lambda x: is_valid_file(parser, x),
    )
    parser.add_argument(
        "--copy",
        action="store_true",
        required="--out-dir" in " ".join(sys.argv),
    )
    parser.add_argument(
        "-d", "--debug", action="store_true", help="print debug messages to stderr"
    )
    parser.add_argument(
        "--out-dir",
        dest="out_dir",
        type=str,
        required="--copy" in " ".join(sys.argv),
        help="Directory where dereplicated assemblies will be copied",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help="print debug messages to stderr",
    )
    args = parser.parse_args()

    return args


def create_jobs_db(db, path):
    pass


def load_classifications(tax_file):
    classifications = collections.defaultdict(list)
    for line in tax_file:
        parts = line.strip().split("\t")
        accession = parts[0]
        if accession.startswith("RS_") or accession.startswith("GB_"):
            accession = accession[3:]
        taxon = parts[1]
        classifications[taxon].append(accession)
    return classifications


def find_all_assemblies(in_dir):
    log.info("Finding assemblies in {}".format(in_dir))
    all_assemblies = [
        str(x)
        for x in sorted(pathlib.Path(in_dir).absolute().glob("**/*"))
        if x.is_file()
    ]
    log.info("Found {:,} files in {}".format(len(all_assemblies), in_dir))
    return all_assemblies


def find_assemblies_for_accessions(accessions, all_assemblies):
    acc_to_assemblies = {}
    found_count, total_count = 0, 0
    not_found = []
    for accession in accessions:
        total_count += 1
        assembly_filename = get_assembly_filename(accession, all_assemblies)
        if assembly_filename is not None:
            found_count += 1
            acc_to_assemblies[accession] = assembly_filename
        else:
            not_found.append(accession)

    # logging.info("Found {}/{} assemblies".format(found_count, total_count))

    # if not_found:
    #     logging.info("Failed to find assemblies for the following accessions:")
    #     wrapper = textwrap.TextWrapper(
    #         initial_indent="    ", subsequent_indent="    ", width=100
    #     )
    #     logging.info(wrapper.fill(", ".join(not_found)))

    return acc_to_assemblies


def get_assembly_filename(accession, all_assemblies):
    if accession.startswith("GCF_") or accession.startswith("GCA_"):
        accession = accession.split(".")[0]
        assert len(accession) == 13
    accession_dot = accession + "."
    accession_under = accession + "_"
    assembly_filenames = [
        x
        for x in all_assemblies
        if x.rpartition("/")[2].startswith(accession_dot)
        or x.rpartition("/")[2].startswith(accession_under)
    ]
    if len(assembly_filenames) == 0:
        return None
    elif len(assembly_filenames) > 1:
        sys.exit(
            "\nError: ambiguous assembly filenames, accession={}, "
            "filenames={}".format(accession, assembly_filenames)
        )
    return assembly_filenames[0]


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {
        "gz": (b"\x1f", b"\x8b", b"\x08"),
        "bz2": (b"\x42", b"\x5a", b"\x68"),
        "zip": (b"\x50", b"\x4b", b"\x03", b"\x04"),
    }
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, "rb")
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = "plain"
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == "bz2":
        sys.exit("Error: cannot use bzip2 format - use gzip instead")
        sys.exit("Error: cannot use zip format - use gzip instead")
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == "gz":
        return gzip.open
    else:  # plain text
        return open


def get_assembly_length(filename):
    contig_lengths = sorted(get_contig_lengths(filename), reverse=True)
    total_length = sum(contig_lengths)
    return total_length


def get_contig_lengths(filename):
    lengths = []
    with get_open_func(filename)(filename, "rt") as fasta_file:
        name = ""
        sequence = ""
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":  # Header line = start of new contig
                if name:
                    lengths.append(len(sequence))
                    sequence = ""
                name = line[1:].split()[0]
            else:
                sequence += line
        if name:
            lengths.append(len(sequence))
    return lengths


@contextmanager
def suppress_stdout():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(devnull, "w") as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


def applyParallel(dfGrouped, func, threads, parms):
    p = Pool(threads)
    func = partial(func, parms=parms)
    ret_list = tqdm.tqdm(
        p.map(func, [group for name, group in dfGrouped]),
        total=len([group for name, group in dfGrouped]),
    )
    return pd.concat(ret_list)