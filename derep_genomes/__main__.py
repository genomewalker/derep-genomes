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
import tempfile
import os
from derep_genomes.general import (
    get_arguments,
    find_all_assemblies,
    load_classifications,
    process_one_taxon,
)


def main():
    args = get_arguments()
    all_assemblies = find_all_assemblies(args.in_dir)
    os.makedirs(args.out_dir, exist_ok=True)
    classifications = load_classifications(args.tax_file)
    for taxon in sorted(classifications.keys()):
        process_one_taxon(
            taxon,
            classifications[taxon],
            all_assemblies,
            args.out_dir,
            args.threads,
            args.threshold,
        )
    print()


if __name__ == "__main__":
    main()
