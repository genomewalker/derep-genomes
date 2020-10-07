"""
This script will dereplicate GTDB assemblies to a user-specified threshold, copying the
dereplicated assemblies to a new directory.

Usage:
    dereplicate_assemblies.py --threshold 0.005 assemblies derep bac_and_arc_taxonomy_r86.tsv

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

from derep_genomes import __version__

import subprocess
import sys
import os
import shutil
from derep_genomes.general import (
    get_arguments,
    find_all_assemblies,
    load_classifications,
    find_assemblies_for_accessions,
)
from derep_genomes.lgraph import dereplicate
import logging
import pathlib

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


def process_one_taxon(
    classification,
    accessions,
    all_assemblies,
    out_dir,
    threads,
    threshold,
    chunks,
    slurm_config,
    tmp_dir,
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
        derep_assemblies = dereplicate(
            acc_to_assemblies, threads, threshold, chunks, slurm_config, tmp_dir
        )
        print("Copying dereplicated assemblies to output directory:")
        for assembly in derep_assemblies:
            print("    {} -> {}".format(assembly, out_dir))
            shutil.copy(assembly, out_dir)


def main():
    args = get_arguments()
    all_assemblies = find_all_assemblies(args.in_dir)
    os.makedirs(args.out_dir, exist_ok=True)
    classifications = load_classifications(args.tax_file)
    tmp_dir = pathlib.Path(args.tmp_dir).absolute()
    for taxon in sorted(classifications.keys()):
        process_one_taxon(
            taxon,
            classifications[taxon],
            all_assemblies,
            args.out_dir,
            args.threads,
            args.threshold,
            args.chunks,
            args.slurm_config,
            tmp_dir,
        )


if __name__ == "__main__":
    main()
