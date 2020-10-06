import argparse
import sys
import shutil
import collections
import gzip
import os
import pathlib
import collections
import textwrap
from derep_genomes.lgraph import dereplicate


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


def process_one_taxon(
    classification, accessions, all_assemblies, out_dir, threads, threshold
):
    accessions = sorted(accessions)
    print(classification)
    acc_to_assemblies = find_assemblies_for_accessions(accessions, all_assemblies)
    if len(acc_to_assemblies) == 0:
        return
    if len(acc_to_assemblies) == 1:
        only_assembly = list(acc_to_assemblies.values())[0]
        print("Only one assembly for this species, copying to output directory:")
        print("    {} -> {}".format(only_assembly, out_dir))
        shutil.copy(only_assembly, out_dir)
    else:
        print(
            "{:,} assemblies for this species, clustering to dereplicate.".format(
                len(acc_to_assemblies)
            )
        )
        derep_assemblies = dereplicate(acc_to_assemblies, threads, threshold)
        print("Copying dereplicated assemblies to output directory:")
        for assembly in derep_assemblies:
            print("    {} -> {}".format(assembly, out_dir))
            shutil.copy(assembly, out_dir)


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
        "--threshold",
        type=float,
        default=2,
        help="Z-score filtering threshold",
    )
    parser.add_argument(
        "--threads", type=int, default=16, help="Number of threads (for Mash)"
    )
    args = parser.parse_args()
    return args
