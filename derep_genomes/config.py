"""
Configuration management for genome dereplication.
"""

import os
import sys
import time
import yaml
import argparse
from pathlib import Path
from typing import Dict, Any, Optional, List, Union

import derep_genomes
from derep_genomes.utils.helpers import check_command_exists, safe_file_path


class Config:
    """Configuration class for genome dereplication."""

    DEFAULT_VALUES = {
        # General parameters
        "db_path": "derep-genomes.db",
        "threads": 2,
        "workers": 2,
        "slurm_threads": 2,
        "tmp_dir": "/tmp",
        "threshold": 2.0,
        "batch_size": 1000,
        "log_level": "info",
        "log_file": None,
        "similarity_method": "fastani",
        "selected_taxa": None,
        "tax_file": None,
        "prefix": time.strftime("%Y%m%d-%H%M%S"),
        "out_dir": None,
        "slurm_config": None,
        "assm_max": 10,
        "chunks": 10,
        "max_jobs_array": 1000,
        # FastANI parameters
        "ani_fraglen": None,  # Auto-calculated when None
        "ani_fraglen_fraction": 0.005,
        "ani_min_frac": 0.0,
        "ani_threshold": 95.0,
        "min_genome_size": 20000000,
        # SKANI parameters
        "skani_k": 15,
        "skani_s": 10000,
        "skani_min_conf": 0.2,
        "skani_target_ani": 0.95,
        "skani_min_kmers": 100,
        # Mash parameters
        "mash_s": 10000,
        "mash_k": 21,
        "mash_threshold": 0.01,
        # Dashing parameters
        "dashing_k": 21,
        "dashing_s": 12,
        "dashing_threshold": 0.01,
        # Common pre-clustering parameters (for both Mash and Dashing)
        "xash_threshold": 0.01,
    }

    def __init__(self):
        """Initialize configuration with default values."""
        # Set default values
        for key, value in self.DEFAULT_VALUES.items():
            setattr(self, key, value)

        # Generate timestamp for file names
        self.timestamp = time.strftime("%Y%m%d-%H%M%S")
        self.prefix = self.timestamp

        # Flags
        self.debug = False
        self.copy = False
        self.dashing = False
        self.skani = False

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "Config":
        """
        Create configuration from command-line arguments.

        Args:
            args: Command-line arguments

        Returns:
            Configuration object
        """
        config = cls()

        # Copy all attributes from args to config
        for key, value in vars(args).items():
            if value is not None:  # Only set non-None values
                setattr(config, key, value)

        # Convert paths to absolute paths
        if config.tmp_dir:
            config.tmp_dir = safe_file_path(config.tmp_dir)

        if config.out_dir:
            config.out_dir = safe_file_path(config.out_dir)

        if config.selected_taxa:
            config.selected_taxa = safe_file_path(config.selected_taxa)

        if config.tax_file:
            config.tax_file = safe_file_path(config.tax_file)

        if config.slurm_config:
            config.slurm_config = safe_file_path(config.slurm_config)

        # Set log level
        if config.debug:
            config.log_level = "debug"

        # Set similarity method
        if config.skani:
            config.similarity_method = "skani"

        # Check for required commands
        config.check_required_commands()

        return config

    @classmethod
    def from_file(cls, file_path: str) -> "Config":
        """
        Load configuration from a YAML file.

        Args:
            file_path: Path to configuration file

        Returns:
            Configuration object
        """
        config = cls()

        try:
            with open(file_path, "r") as f:
                config_dict = yaml.safe_load(f)

            # Update configuration
            for key, value in config_dict.items():
                if hasattr(config, key):
                    setattr(config, key, value)

            # Convert paths to absolute paths
            if config.tmp_dir:
                config.tmp_dir = safe_file_path(config.tmp_dir)

            if config.out_dir:
                config.out_dir = safe_file_path(config.out_dir)

            if config.selected_taxa:
                config.selected_taxa = safe_file_path(config.selected_taxa)

            if config.tax_file:
                config.tax_file = safe_file_path(config.tax_file)

            if config.slurm_config:
                config.slurm_config = safe_file_path(config.slurm_config)

            # Set log level
            if config.debug:
                config.log_level = "debug"

            # Set similarity method
            if config.skani:
                config.similarity_method = "skani"

            # Check for required commands
            config.check_required_commands()

        except Exception as e:
            print(f"Error loading configuration from {file_path}: {e}")
            sys.exit(1)

        return config

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert configuration to dictionary.

        Returns:
            Dictionary representation of configuration
        """
        return {
            key: getattr(self, key)
            for key in dir(self)
            if not key.startswith("_") and not callable(getattr(self, key))
        }

    def save(self, file_path: str) -> None:
        """
        Save configuration to a YAML file.

        Args:
            file_path: Path to output file
        """
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(file_path), exist_ok=True)

        # Save configuration
        with open(file_path, "w") as f:
            yaml.dump(self.to_dict(), f)

    def check_required_commands(self) -> None:
        """Check if required commands are available."""
        required_commands = ["fastANI", "mash"]

        if self.dashing:
            required_commands.append("dashing")

        if self.skani:
            required_commands.append("skani")

        for command in required_commands:
            if not check_command_exists(command):
                print(f"Error: {command} not found in PATH")
                sys.exit(1)

    def __str__(self) -> str:
        """String representation of configuration."""
        lines = [f"Configuration:"]
        for key, value in sorted(self.to_dict().items()):
            if key != "DEFAULT_VALUES":
                lines.append(f"  {key}: {value}")
        return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Dereplicate genomes based on similarity",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    timestamp = time.strftime("%Y%m%d-%H%M%S")

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--data",
        required=True,
        metavar="FILE",
        dest="tax_file",
        help="TSV file with the genome information",
    )

    # Optional arguments
    parser.add_argument(
        "--db",
        type=str,
        metavar="DB",
        default="derep-genomes.db",
        help="Path to DuckDB database file",
    )
    parser.add_argument(
        "--workers",
        type=int,
        metavar="INT",
        default=2,
        help="Number of worker processes",
    )
    parser.add_argument(
        "--threads",
        type=int,
        metavar="INT",
        default=2,
        help="Number of threads per worker",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default=timestamp,
        metavar="PREFIX",
        help="Prefix for the file name results",
    )
    parser.add_argument(
        "--tmp",
        type=str,
        default="/tmp",
        metavar="DIR",
        dest="tmp_dir",
        help="Temporary directory",
    )
    parser.add_argument(
        "--threshold",
        metavar="FLOAT",
        type=float,
        default=2.0,
        help="Z-score filtering threshold",
    )
    parser.add_argument(
        "--max-assemblies",
        type=int,
        metavar="INT",
        dest="assm_max",
        default=10,
        help="Maximum number of assemblies to process for small batches",
    )
    parser.add_argument(
        "--selected-taxa",
        dest="selected_taxa",
        metavar="FILE",
        help="File with selected taxa. Taxa should have the format: d__;p__;c__;o__;f__;g__;s__",
    )
    parser.add_argument(
        "--copy", action="store_true", help="Copy assembly files to the output folder"
    )
    parser.add_argument(
        "--out-dir",
        dest="out_dir",
        type=str,
        help="Directory where dereplicated assemblies will be copied",
    )
    parser.add_argument(
        "--config", metavar="FILE", help="Path to configuration file (YAML)"
    )

    # Method selection groups
    precluster_group = parser.add_argument_group("pre-clustering method")
    precluster_group.add_argument(
        "--dashing",
        action="store_true",
        help="Use Dashing instead of Mash for initial pre-clustering",
    )

    analysis_group = parser.add_argument_group("detailed analysis method")
    analysis_group.add_argument(
        "--skani",
        action="store_true",
        help="Use SKANI instead of FastANI for detailed similarity calculation",
    )

    # FastANI parameters
    fastani_group = parser.add_argument_group("FastANI parameters")
    fastani_group.add_argument(
        "--ani-fraglen",
        type=int,
        metavar="INT",
        help="Fragment length for FastANI (default: auto-calculated)",
    )
    fastani_group.add_argument(
        "--ani-fraglen-fraction",
        type=float,
        metavar="FLOAT",
        default=0.005,
        help="Fraction of genome size for ANI fragment length",
    )
    fastani_group.add_argument(
        "--ani-min-frac",
        type=float,
        metavar="FLOAT",
        default=0.0,
        help="Minimum fraction of genome that must align",
    )
    fastani_group.add_argument(
        "--ani-threshold",
        type=float,
        metavar="FLOAT",
        default=95.0,
        help="Minimum ANI value to consider genomes as similar",
    )
    fastani_group.add_argument(
        "--min-genome-size",
        dest="min_genome_size",
        metavar="INT",
        type=int,
        default=20000000,
        help="Minimum genome size for ANI fragment length heuristics",
    )

    # SKANI parameters
    skani_group = parser.add_argument_group("SKANI parameters")
    skani_group.add_argument(
        "--skani-k",
        type=int,
        metavar="INT",
        default=15,
        help="K-mer size for SKANI sketch creation",
    )
    skani_group.add_argument(
        "--skani-s",
        type=int,
        metavar="INT",
        default=10000,
        help="Sketch size for SKANI",
    )
    skani_group.add_argument(
        "--skani-min-conf",
        type=float,
        metavar="FLOAT",
        default=0.2,
        help="Minimum confidence value for SKANI reporting",
    )
    skani_group.add_argument(
        "--skani-target-ani",
        type=float,
        metavar="FLOAT",
        default=0.95,
        help="Target ANI value for SKANI optimization",
    )
    skani_group.add_argument(
        "--skani-min-kmers",
        type=int,
        metavar="INT",
        default=100,
        help="Minimum number of k-mers two genomes must share in SKANI",
    )

    # Mash parameters
    mash_group = parser.add_argument_group("Mash parameters")
    mash_group.add_argument(
        "--mash-s", type=int, metavar="INT", default=10000, help="Sketch size for Mash"
    )
    mash_group.add_argument(
        "--mash-k",
        type=int,
        metavar="INT",
        default=21,
        help="K-mer size for Mash sketch creation",
    )
    mash_group.add_argument(
        "--mash-threshold",
        type=float,
        metavar="FLOAT",
        default=0.01,
        help="Maximum Mash distance to consider genomes as similar",
    )

    # Dashing parameters
    dashing_group = parser.add_argument_group("Dashing parameters")
    dashing_group.add_argument(
        "--dashing-k",
        type=int,
        metavar="INT",
        default=21,
        help="K-mer size for Dashing sketch creation",
    )
    dashing_group.add_argument(
        "--dashing-s",
        type=int,
        metavar="INT",
        default=12,
        help="Sketch size factor for Dashing",
    )
    dashing_group.add_argument(
        "--dashing-threshold",
        type=float,
        metavar="FLOAT",
        default=0.01,
        help="Maximum Dashing distance to consider genomes as similar",
    )

    # Common pre-clustering parameters
    parser.add_argument(
        "--xash-threshold",
        metavar="FLOAT",
        type=float,
        default=0.01,
        help="Threshold for Mash/Dashing pre-clustering distance",
    )

    # SLURM arguments
    slurm = parser.add_argument_group("SLURM arguments")
    slurm.add_argument(
        "--slurm-config",
        dest="slurm_config",
        help="YAML configuration file for SLURM",
        metavar="FILE",
    )
    slurm.add_argument(
        "--slurm-threads",
        type=int,
        metavar="INT",
        dest="slurm_threads",
        default=2,
        help="Number of threads to use in a SLURM job",
    )
    slurm.add_argument(
        "--chunk-size",
        metavar="INT",
        dest="chunks",
        type=int,
        default=10,
        help="Number of genomes in each chunk for FastANI in SLURM",
    )
    slurm.add_argument(
        "--slurm-arr-size",
        dest="max_jobs_array",
        metavar="INT",
        type=int,
        default=1000,
        help="Slurm maximum job array size",
    )

    # Debug flag
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")

    # Version
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {derep_genomes.__version__}",
        help="Show program version",
    )

    args = parser.parse_args()

    # Validate arguments
    if args.copy and not args.out_dir:
        parser.error("--out-dir is required when --copy is specified")

    if args.out_dir and not args.copy:
        parser.error("--copy is required when --out-dir is specified")

    # Load configuration file if specified
    if args.config:
        # Load config file and override with command-line arguments
        config = Config.from_file(args.config)

        # Override config with command-line arguments
        for key, value in vars(args).items():
            if value is not None and key != "config":
                setattr(config, key, value)

        # Convert config back to args
        for key, value in config.to_dict().items():
            if hasattr(args, key):
                setattr(args, key, value)

    return args
