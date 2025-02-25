"""
Utility functions for genome dereplication.
"""

import os
import gzip
import logging
import subprocess
from pathlib import Path
from typing import List, Dict, Tuple, Set, Union, Optional, Any, Iterator, Callable
from concurrent.futures import ProcessPoolExecutor
import functools

import numpy as np
from Bio import SeqIO

logger = logging.getLogger("derep_genomes")


class GenomeLengthCache:
    """Cache for storing genome lengths to avoid repeated calculations."""

    def __init__(self):
        self.cache = {}
        self.cache_hits = 0
        self.cache_misses = 0

    def get(self, genome_path: str) -> int:
        """
        Get the length of a genome, using cached value if available.

        Args:
            genome_path: Path to genome file

        Returns:
            Genome length in base pairs
        """
        if genome_path in self.cache:
            self.cache_hits += 1
            return self.cache[genome_path]

        # Calculate length
        self.cache_misses += 1
        length = self._calculate_length(genome_path)
        self.cache[genome_path] = length
        return length

    def _calculate_length(self, genome_path: str) -> int:
        """
        Calculate the length of a genome.

        Args:
            genome_path: Path to genome file

        Returns:
            Genome length in base pairs
        """
        total_length = 0
        try:
            with get_open_func(genome_path)(genome_path, "rt") as fasta_file:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    total_length += len(record.seq)
            return total_length
        except Exception as e:
            logger.warning(f"Error getting assembly length for {genome_path}: {e}")
            return 0

    def prefetch(self, genome_paths: List[str], max_workers: int = 16) -> None:
        """
        Prefetch genome lengths in parallel.

        Args:
            genome_paths: List of genome file paths
            max_workers: Maximum number of parallel workers
        """
        missing_paths = [p for p in genome_paths if p not in self.cache]
        if not missing_paths:
            return

        # Use parallel processing to calculate lengths
        with ProcessPoolExecutor(
            max_workers=min(max_workers, len(missing_paths))
        ) as executor:
            lengths = list(executor.map(self._calculate_length, missing_paths))

        # Update cache
        for path, length in zip(missing_paths, lengths):
            self.cache[path] = length

    def get_stats(self) -> Dict[str, Any]:
        """
        Get cache statistics.

        Returns:
            Dictionary with cache statistics
        """
        return {
            "hits": self.cache_hits,
            "misses": self.cache_misses,
            "size": len(self.cache),
            "hit_ratio": (
                self.cache_hits / (self.cache_hits + self.cache_misses)
                if (self.cache_hits + self.cache_misses) > 0
                else 0
            ),
        }


# Initialize cache at module level
genome_length_cache = GenomeLengthCache()


def get_compression_type(filename: str) -> str:
    """
    Determine the compression type of a file.

    Args:
        filename: Path to the file

    Returns:
        Compression type ('gz', 'bz2', 'zip', or 'plain')
    """
    magic_dict = {
        "gz": (b"\x1f", b"\x8b", b"\x08"),
        "bz2": (b"\x42", b"\x5a", b"\x68"),
        "zip": (b"\x50", b"\x4b", b"\x03", b"\x04"),
    }

    max_len = max(len(x) for x in magic_dict)

    try:
        with open(filename, "rb") as f:
            file_start = f.read(max_len)

        for file_type, magic_bytes in magic_dict.items():
            if file_start.startswith(magic_bytes):
                return file_type
    except Exception as e:
        logger.warning(f"Could not determine compression type of {filename}: {e}")

    return "plain"


def get_open_func(filename: str) -> Callable:
    """
    Get the appropriate function to open a file based on its compression.

    Args:
        filename: Path to the file

    Returns:
        Function to open the file
    """
    compression = get_compression_type(filename)

    if compression == "gz":
        return gzip.open
    else:  # plain text or other
        return open


def get_assembly_length(filename: str) -> int:
    """
    Get the total length of an assembly (cached).

    Args:
        filename: Path to the assembly file

    Returns:
        Total length of the assembly
    """
    return genome_length_cache.get(filename)


def get_contig_lengths(filename: str) -> List[int]:
    """
    Get the lengths of all contigs in an assembly.

    Args:
        filename: Path to the assembly file

    Returns:
        List of contig lengths
    """
    lengths = []

    try:
        with get_open_func(filename)(filename, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                lengths.append(len(record.seq))
    except Exception as e:
        logger.warning(f"Error getting contig lengths for {filename}: {e}")
        return []

    return lengths


def get_assembly_n50(filename: str) -> int:
    """
    Get the N50 of an assembly.

    Args:
        filename: Path to the assembly file

    Returns:
        N50 of the assembly
    """
    contig_lengths = sorted(get_contig_lengths(filename), reverse=True)

    if not contig_lengths:
        return 0

    total_length = sum(contig_lengths)
    target_length = total_length * 0.5
    length_so_far = 0

    for contig_length in contig_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            return contig_length

    return 0


def check_command_exists(command: str) -> bool:
    """
    Check if a command exists in the system PATH.

    Args:
        command: Command to check

    Returns:
        True if the command exists, False otherwise
    """
    try:
        subprocess.run(
            [command, "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        return True
    except (subprocess.SubprocessError, FileNotFoundError):
        return False


def get_accession_from_path(path: str) -> str:
    """
    Extract the accession from a file path.

    Args:
        path: File path

    Returns:
        Accession
    """
    basename = os.path.basename(path)

    # Extract accession based on common formats
    if basename.startswith("GCF_") or basename.startswith("GCA_"):
        accession = basename.split(".")[0]
        if len(accession) == 13:  # Standard NCBI accession length
            return accession

    # Try other patterns
    parts = basename.split("_")
    if len(parts) > 1:
        return parts[0]

    # Default case, remove extension
    return os.path.splitext(basename)[0]


def find_assembly_file(accession: str, directory: str) -> Optional[str]:
    """
    Find the assembly file for an accession in a directory.

    Args:
        accession: Accession to find
        directory: Directory to search

    Returns:
        Path to the assembly file, or None if not found
    """
    # Normalize accession
    if accession.startswith("GCF_") or accession.startswith("GCA_"):
        accession = accession.split(".")[0]

    # Patterns to search for
    patterns = [f"{accession}.", f"{accession}_"]

    # Search for files matching patterns
    for root, _, files in os.walk(directory):
        for filename in files:
            for pattern in patterns:
                if filename.startswith(pattern):
                    return os.path.join(root, filename)

    return None


def split_into_chunks(items: List[Any], chunk_size: int) -> List[List[Any]]:
    """
    Split a list into chunks of fixed size.

    Args:
        items: List to split
        chunk_size: Maximum chunk size

    Returns:
        List of chunks
    """
    return [items[i : i + chunk_size] for i in range(0, len(items), chunk_size)]


def safe_file_path(path: str) -> str:
    """
    Ensure a file path is safe and absolute.

    Args:
        path: File path

    Returns:
        Safe absolute file path
    """
    return str(Path(path).absolute())


def sanitize_for_filename(text: str) -> str:
    """
    Sanitize a string to make it safe for use in a filename.

    Args:
        text: Input string

    Returns:
        Sanitized string
    """
    # Replace problematic characters
    replacements = {
        ";": "_",
        " ": "_",
        ":": "_",
        "/": "_",
        "\\": "_",
        "*": "_",
        "?": "_",
        '"': "_",
        "<": "_",
        ">": "_",
        "|": "_",
    }

    for char, replacement in replacements.items():
        text = text.replace(char, replacement)

    return text


def count_file_lines(file_path: str) -> int:
    """
    Count the number of lines in a file efficiently.

    Args:
        file_path: Path to the file

    Returns:
        Number of lines
    """
    with open(file_path, "rb") as f:
        return sum(1 for _ in f)
