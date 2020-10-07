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

from derep_genomes.__init__ import __version__

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


def get_arguments():
    parser = argparse.ArgumentParser(description="Cluster assemblies in each taxon")
    parser.add_argument(
        "in_dir", type=str, help="Directory containing all GTDB assemblies"
    )
    parser.add_argument(
        "out_dir",
        type=str,
        help="Directory where dereplicated assemblies will be copied",
    )
    parser.add_argument("tax_file", type=str, help="GTDB taxonomy file")
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
        "--chunks",
        type=int,
        default=32,
        help="Number of chunks to create (for fastANI in slurm)",
    )
    parser.add_argument(
        "-d", "--debug", action="store_true", help="print debug messages to stderr"
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


def load_classifications(tax_file):
    classifications = collections.defaultdict(list)
    with open(tax_file, "rt") as tax:
        for line in tax:
            parts = line.strip().split("\t")
            accession = parts[0]
            if accession.startswith("RS_") or accession.startswith("GB_"):
                accession = accession[3:]
            taxon = parts[1]
            classifications[taxon].append(accession)
    return classifications


def find_all_assemblies(in_dir):
    print("\nLooking for files in {}:".format(in_dir))
    all_assemblies = [
        str(x) for x in sorted(pathlib.Path(in_dir).glob("**/*")) if x.is_file()
    ]
    print("found {:,} files".format(len(all_assemblies)))
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
        print(
            "\r{:,} / {:,} assemblies found".format(found_count, total_count),
            end="",
            flush=True,
        )
    if not_found:
        print("    failed to find assemblies for the following accessions:")
        wrapper = textwrap.TextWrapper(
            initial_indent="    ", subsequent_indent="    ", width=100
        )
        print(wrapper.fill(", ".join(not_found)))
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