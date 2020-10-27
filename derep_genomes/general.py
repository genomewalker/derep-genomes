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


def check_values(val, minval, maxval, parser, var):
    value = float(val)
    if value < minval or value > maxval:
        parser.error(
            "argument %s: Invalid value %s. Range has to be between %s and %s!"
            % (
                var,
                value,
                minval,
                maxval,
            )
        )
    return value


# From: https://stackoverflow.com/a/11541450
def is_valid_file(parser, arg, var):
    if not os.path.exists(arg):
        if os.path.isfile(arg):
            parser.error("argument %s: The file %s does not exist!" % (var, arg))
        else:
            parser.error("argument %s: The directory %s does not exist!" % (var, arg))
    else:
        if os.path.isfile(arg):
            return open(arg, "r")  # return an open file handle
        else:
            return arg


help_msg = {
    "in_dir": "Directory containing all assemblies",
    "tax_file": "TSV file with the taxonomic information",
    "db": "SQLite3 DB to store the results",
    "tmp": "Temporary directory",
    "threads": "Number of threads (for fastANI)",
    "threshold": "Z-score filtering threshold",
    "slurm_config": "YAML configuration file for SLURM",
    "chunks": "Number of genomes in each chunk for fastANI in SLURM",
    "max_jobs_array": "Slurm maximum job array size",
    "selected_taxa": "File with selected taxa. Taxa should have the following formar: d__;p__;c__;o__;f__;g__;s__",
    "out_dir": "Directory where dereplicated assemblies will be copied",
    "debug": "Print debug messages",
    "version": "Print program version",
    "mash_threshold": "Mash distance threshold where to filter",
}


def get_arguments(argv=None):

    parser = argparse.ArgumentParser(
        description="Cluster assemblies in each taxon",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    optional = parser._action_groups.pop()  # Edited this line
    required = parser.add_argument_group("required arguments")
    parser._action_groups.append(optional)  # added this line
    slurm = parser.add_argument_group("SLURM arguments")

    required.add_argument(
        "--in-dir",
        dest="in_dir",
        default=argparse.SUPPRESS,
        required=True,
        metavar="DIR",
        type=lambda x: is_valid_file(parser, x, "--in-dir"),
        help=help_msg["in_dir"],
    )
    required.add_argument(
        "--taxa",
        required=True,
        metavar="FILE",
        default=argparse.SUPPRESS,
        type=lambda x: is_valid_file(parser, x, "--taxa"),
        dest="tax_file",
        help=help_msg["tax_file"],
    )
    optional.add_argument(
        "--db",
        dest="db",
        type=str,
        metavar="DB",
        default="derep-genomes.db",
        help=help_msg["db"],
    )
    optional.add_argument(
        "--threads",
        type=int,
        metavar="INT",
        dest="threads",
        default=16,
        help=help_msg["threads"],
    )
    optional.add_argument(
        "--tmp",
        type=str,
        default="/tmp",
        metavar="DIR",
        dest="tmp_dir",
        help=help_msg["tmp"],
    )
    optional.add_argument(
        "--threshold",
        metavar="FLOAT",
        type=float,
        default=2.0,
        help=help_msg["threshold"],
    )
    optional.add_argument(
        "--mash-threshold",
        metavar="FLOAT",
        type=lambda x: check_values(
            x, minval=0, maxval=1, parser=parser, var="--mash-threshold"
        ),
        default=0.01,
        help=help_msg["mash_threshold"],
        dest="mash_threshold",
    )
    slurm.add_argument(
        "--slurm-config",
        dest="slurm_config",
        default=argparse.SUPPRESS,
        required=False,
        help=help_msg["slurm_config"],
        metavar="FILE",
        type=lambda x: is_valid_file(parser, x, "--slurm-config"),
    )
    slurm.add_argument(
        "--chunk-size",
        metavar="INT",
        dest="chunks",
        type=int,
        default=10,
        help=help_msg["chunks"],
    )
    slurm.add_argument(
        "--slurm-arr-size",
        dest="max_jobs_array",
        metavar="INT",
        type=int,
        default=1000,
        help=help_msg["max_jobs_array"],
    )
    optional.add_argument(
        "--selected-taxa",
        dest="selected_taxa",
        required=False,
        metavar="FILE",
        default=argparse.SUPPRESS,
        type=lambda x: is_valid_file(parser, x, "--selected-taxa"),
        help=help_msg["selected_taxa"],
    )
    optional.add_argument(
        "--copy",
        action="store_true",
        required="--out-dir" in " ".join(sys.argv),
        help="Copy assembly files to the output folder",
    )
    optional.add_argument(
        "--out-dir",
        dest="out_dir",
        default=argparse.SUPPRESS,
        type=str,
        required="--copy" in " ".join(sys.argv),
        help=help_msg["out_dir"],
    )
    optional.add_argument(
        "--debug", dest="debug", action="store_true", help=help_msg["debug"]
    )
    optional.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help=help_msg["version"],
    )
    args = parser.parse_args()
    if not hasattr(args, "out_dir"):
        args.out_dir = None
    if not hasattr(args, "selected_taxa"):
        args.selected_taxa = None
    if not hasattr(args, "slurm_config"):
        args.slurm_config = None

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


# from https://stackoverflow.com/a/62478211
def absolute_file_paths(directory):
    path = os.path.abspath(directory)
    files = list(os.scandir(path))
    return [entry.path for entry in files if entry.is_file()]


def find_all_assemblies(in_dir):
    log.info("Finding assemblies in {}".format(in_dir))
    all_assemblies = absolute_file_paths(in_dir)
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


def get_assembly_n50(filename):
    contig_lengths = sorted(get_contig_lengths(filename), reverse=True)
    total_length = sum(contig_lengths)
    target_length = total_length * 0.5
    length_so_far = 0
    for contig_length in contig_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            return contig_length
    return 0


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