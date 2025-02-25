"""
Base class for genome similarity measures.
"""

import abc
import os
import tempfile
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Set, Union, Optional, Any

import numpy as np

from derep_genomes.utils.helpers import genome_length_cache, sanitize_for_filename

logger = logging.getLogger(__name__)


class SimilarityMeasure(abc.ABC):
    """Abstract base class for genome similarity measures."""

    def __init__(self, threads: int, threshold: float, temp_dir: Path):
        """
        Initialize the similarity measure.

        Args:
            threads: Number of threads to use
            threshold: Threshold for similarity
            temp_dir: Directory for temporary files
        """
        self.threads = threads
        self.threshold = threshold
        self.temp_dir = temp_dir

        # Ensure temp directory exists
        os.makedirs(self.temp_dir, exist_ok=True)

    @abc.abstractmethod
    def compute_pairwise(self, assemblies: List[str]) -> Dict[str, Any]:
        """
        Compute pairwise similarities between assemblies.

        Args:
            assemblies: List of assembly file paths

        Returns:
            Dictionary with pairwise similarity data
        """
        pass

    @abc.abstractmethod
    def filter_by_similarity(self, similarity_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Filter pairwise similarities by threshold.

        Args:
            similarity_data: Dictionary with pairwise similarity data

        Returns:
            Filtered similarity data
        """
        pass

    @abc.abstractmethod
    def create_similarity_graph(self, similarity_data: Dict[str, Any]) -> Any:
        """
        Create a similarity graph from pairwise similarities.

        Args:
            similarity_data: Dictionary with pairwise similarity data

        Returns:
            Graph representation of similarities
        """
        pass

    def create_temp_file(
        self, suffix: str = ""
    ) -> Tuple[str, tempfile.NamedTemporaryFile]:
        """
        Create a temporary file in the temp directory.

        Args:
            suffix: File suffix

        Returns:
            Tuple of (file path, file object)
        """
        temp_file = tempfile.NamedTemporaryFile(
            dir=self.temp_dir, prefix="derep-", suffix=suffix, delete=False, mode="w+t"
        )
        return temp_file.name, temp_file

    def write_assembly_list(self, assemblies: List[str]) -> str:
        """
        Write assembly list to a temporary file.

        Args:
            assemblies: List of assembly file paths

        Returns:
            Path to the temporary file
        """
        path, tmp_file = self.create_temp_file(suffix=".txt")
        try:
            for assembly in assemblies:
                tmp_file.write(f"{assembly}\n")
            tmp_file.flush()
        finally:
            tmp_file.close()
        return path

    def check_missing_assemblies(
        self, similarity_data: Dict[str, Any], all_assemblies: List[str]
    ) -> List[str]:
        """
        Check which assemblies are missing from the similarity data.

        Args:
            similarity_data: Dictionary with pairwise similarity data
            all_assemblies: List of all assembly file paths

        Returns:
            List of missing assembly file paths
        """
        # Get unique sources and targets from similarity data
        present_assemblies = set()
        for source, target in zip(
            similarity_data.get("source", []), similarity_data.get("target", [])
        ):
            present_assemblies.add(source)
            present_assemblies.add(target)

        # Find missing assemblies
        missing = [
            assembly
            for assembly in all_assemblies
            if assembly not in present_assemblies
        ]

        if missing:
            logger.debug(
                f"Missing {len(missing)} assemblies in the similarity comparison, "
                "possibly too divergent"
            )

        return missing

    def prefetch_genome_lengths(self, assemblies: List[str]) -> None:
        """
        Prefetch genome lengths for a list of assemblies.

        Args:
            assemblies: List of assembly file paths
        """
        genome_length_cache.prefetch(assemblies)

    def create_taxon_temp_dir(self, taxon: str) -> Path:
        """
        Create a taxon-specific temporary directory.

        Args:
            taxon: Taxonomy name

        Returns:
            Path to the temporary directory
        """
        taxon_safe = sanitize_for_filename(taxon)
        taxon_temp_dir = Path(self.temp_dir) / taxon_safe
        os.makedirs(taxon_temp_dir, exist_ok=True)
        return taxon_temp_dir

    def cleanup_temp_files(self, *files: str) -> None:
        """
        Clean up temporary files.

        Args:
            *files: Paths to files to clean up
        """
        for file_path in files:
            try:
                if file_path and os.path.exists(file_path):
                    os.remove(file_path)
            except OSError as e:
                logger.debug(f"Failed to remove temporary file {file_path}: {e}")
