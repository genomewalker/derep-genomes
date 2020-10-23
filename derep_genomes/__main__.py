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
from derep_genomes.derep import derep
from derep_genomes.index import index

import sys
import argparse
from derep_genomes.general import is_valid_file, help_msg
import logging

log = logging.getLogger("my_logger")


def main():
    logging.basicConfig(
        level=logging.DEBUG, format="%(levelname)s ::: %(asctime)s ::: %(message)s"
    )

    parser = argparse.ArgumentParser(
        description="Cluster assemblies in each taxon",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(help="Functions")
    parser_index = subparsers.add_parser("index", help="...")
    parser_index.set_defaults(parser_index=True, func=index)
    optional_index = parser_index._action_groups.pop()  # Edited this line
    required_index = parser_index.add_argument_group("required arguments")
    required_index.add_argument(
        "--in-dir",
        dest="in_dir",
        default=argparse.SUPPRESS,
        required=True,
        metavar="DIR",
        type=lambda x: is_valid_file(parser_index, x),
        help=help_msg["in_dir"],
    )
    required_index.add_argument(
        "--taxa",
        required=True,
        metavar="FILE",
        default=argparse.SUPPRESS,
        type=lambda x: is_valid_file(parser_index, x),
        dest="tax_file",
        help=help_msg["tax_file"],
    )
    optional_index.add_argument(
        "--threads",
        type=int,
        metavar="INT",
        dest="threads",
        default=16,
        help=help_msg["threads"],
    )
    optional_index.add_argument(
        "--tmp",
        type=str,
        default="/tmp",
        metavar="DIR",
        dest="tmp_dir",
        help=help_msg["tmp"],
    )

    parser_derep = subparsers.add_parser("derep", help="...")
    parser_derep.set_defaults(parser_derep=True, func=derep)
    optional_derep = parser_derep._action_groups.pop()  # Edited this line
    required_derep = parser_derep.add_argument_group("required arguments")
    parser_derep._action_groups.append(optional_derep)  # added this line
    slurm_derep = parser_derep.add_argument_group("SLURM arguments")

    required_derep.add_argument(
        "--in-file",
        dest="in_file",
        default=argparse.SUPPRESS,
        required=True,
        metavar="FILE",
        type=lambda x: is_valid_file(parser_derep, x),
        help=help_msg["in_file"],
    )
    required_derep.add_argument(
        "--taxa",
        required=True,
        metavar="FILE",
        default=argparse.SUPPRESS,
        type=lambda x: is_valid_file(parser_derep, x),
        dest="tax_file",
        help=help_msg["tax_file"],
    )
    optional_derep.add_argument(
        "--db",
        dest="db",
        type=str,
        metavar="DB",
        default="derep-genomes.db",
        help=help_msg["db"],
    )
    optional_derep.add_argument(
        "--threads",
        type=int,
        metavar="INT",
        dest="threads",
        default=16,
        help=help_msg["threads"],
    )
    optional_derep.add_argument(
        "--tmp",
        type=str,
        default="/tmp",
        metavar="DIR",
        dest="tmp_dir",
        help=help_msg["tmp"],
    )
    optional_derep.add_argument(
        "--threshold",
        metavar="FLOAT",
        type=float,
        default=2.0,
        help=help_msg["threshold"],
    )
    slurm_derep.add_argument(
        "--slurm-config",
        dest="slurm_config",
        default=argparse.SUPPRESS,
        required=False,
        help=help_msg["slurm_config"],
        metavar="FILE",
        type=lambda x: is_valid_file(parser_derep, x),
    )
    slurm_derep.add_argument(
        "--chunk-size",
        metavar="INT",
        dest="chunks",
        type=int,
        default=10,
        help=help_msg["chunks"],
    )
    slurm_derep.add_argument(
        "--slurm-arr-size",
        dest="max_jobs_array",
        metavar="INT",
        type=int,
        default=1000,
        help=help_msg["max_jobs_array"],
    )
    optional_derep.add_argument(
        "--selected-taxa",
        dest="selected_taxa",
        required=False,
        metavar="FILE",
        default=argparse.SUPPRESS,
        type=lambda x: is_valid_file(parser_derep, x),
        help=help_msg["selected_taxa"],
    )
    optional_derep.add_argument(
        "--copy",
        action="store_true",
        required="--out-dir" in " ".join(sys.argv),
        help="Copy assembly files to the output folder",
    )
    optional_derep.add_argument(
        "--out-dir",
        dest="out_dir",
        default=argparse.SUPPRESS,
        type=str,
        required="--copy" in " ".join(sys.argv),
        help=help_msg["out_dir"],
    )
    optional_derep.add_argument(
        "--debug", dest="debug", action="store_true", help=help_msg["debug"]
    )
    optional_derep.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help=help_msg["version"],
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    if not hasattr(args, "out_dir"):
        args.out_dir = None
    if not hasattr(args, "selected_taxa"):
        args.selected_taxa = None
    if not hasattr(args, "slurm_config"):
        args.slurm_config = None

    # args = get_arguments()

    args.func(args)


if __name__ == "__main__":
    main()