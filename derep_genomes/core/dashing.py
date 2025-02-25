"""
Dashing implementation for genome similarity pre-clustering.
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Set, Union, Optional, Any

import numpy as np

from derep_genomes.core.similarity import SimilarityMeasure
from derep_genomes.utils.helpers import (
    get_assembly_length,
    genome_length_cache,
    sanitize_for_filename,
)

logger = logging.getLogger(__name__)


class Dashing(SimilarityMeasure):
    """Dashing-based genome similarity measure for pre-clustering."""

    def __init__(
        self,
        threads: int,
        threshold: float,
        temp_dir: Path,
        dashing_k: int = 21,
        dashing_s: int = 12,
        dashing_threshold: float = 0.01,
    ):
        """
        Initialize Dashing similarity measure.

        Args:
            threads: Number of threads to use
            threshold: Threshold for similarity (0-1)
            temp_dir: Directory for temporary files
            dashing_k: K-mer size for Dashing
            dashing_s: Sketch size factor for Dashing
            dashing_threshold: Maximum Dashing distance to consider genomes as similar
        """
        super().__init__(threads, threshold, temp_dir)
        self.dashing_k = dashing_k
        self.dashing_s = dashing_s
        self.dashing_threshold = dashing_threshold

    def run_dashing(self, assembly_list_path: str, output_path: str) -> None:
        """
        Run Dashing on a list of assemblies.

        Args:
            assembly_list_path: Path to file with list of assemblies
            output_path: Path to output file

        Raises:
            subprocess.SubprocessError: If Dashing fails
        """
        cmd = [
            "dashing",
            "dist",
            "-p",
            str(self.threads),
            "-k",
            str(self.dashing_k),
            "-S",
            str(self.dashing_s),
            "--full-mash-dist",
            "--use-bb-minhash",
            "-b",
            "-T",
            "-F",
            assembly_list_path,
            "-M",
            "-O",
            output_path,
        ]

        try:
            logger.debug(f"Running Dashing with command: {' '.join(cmd)}")
            subprocess.run(
                cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )
        except subprocess.SubprocessError as e:
            logger.error(f"Dashing failed: {e}")
            raise

    def compute_pairwise(self, assemblies: List[str]) -> Dict[str, Any]:
        """
        Compute pairwise Dashing similarities between assemblies.

        Args:
            assemblies: List of assembly file paths

        Returns:
            Dictionary with pairwise similarity data
        """
        # Prefetch all genome lengths at once
        genome_length_cache.prefetch(assemblies)

        # Create temporary files
        assembly_list_path = self.write_assembly_list(assemblies)
        output_path, _ = self.create_temp_file(suffix=".dashing.tsv")

        try:
            # Run Dashing
            self.run_dashing(assembly_list_path, output_path)

            # Create graph directly from Dashing output
            from derep_genomes.core.graph import create_graph_from_distance_matrix

            graph = create_graph_from_distance_matrix(
                output_path, self.dashing_threshold
            )

            # Convert graph to edge list for consistent interface
            results = self._graph_to_edge_list(graph)

            return results
        finally:
            # Clean up temporary files
            try:
                os.remove(assembly_list_path)
                if os.path.exists(output_path):
                    os.remove(output_path)
            except OSError:
                pass

    def _graph_to_edge_list(self, graph) -> Dict[str, Any]:
        """
        Convert NetworkX graph to edge list format.

        Args:
            graph: NetworkX graph

        Returns:
            Dictionary with edge list data
        """
        results = {
            "source": [],
            "target": [],
            "distance": [],
            "weight": [],
            "weight_raw": [],
            "source_len": [],
            "target_len": [],
        }

        for u, v, data in graph.edges(data=True):
            results["source"].append(u)
            results["target"].append(v)
            results["weight"].append(data["weight"])
            results["weight_raw"].append(data["weight_raw"])
            results["distance"].append(data.get("distance", 1.0 - data["weight"]))

            results["source_len"].append(graph.nodes[u]["genome_len"])
            results["target_len"].append(graph.nodes[v]["genome_len"])

        return results

    def filter_by_similarity(self, similarity_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Filter Dashing results by distance threshold.

        Args:
            similarity_data: Dictionary with Dashing results

        Returns:
            Filtered similarity data
        """
        # Find indices of results where distance is below threshold
        indices = [
            i
            for i, distance in enumerate(similarity_data.get("distance", []))
            if distance <= self.dashing_threshold
        ]

        # Filter all keys by indices
        filtered_data = {
            key: [similarity_data[key][i] for i in indices] for key in similarity_data
        }

        logger.debug(
            f"Filtered from {len(similarity_data.get('distance', []))} to {len(filtered_data.get('distance', []))} "
            f"comparisons using distance threshold {self.dashing_threshold}"
        )

        return filtered_data

    def create_similarity_graph(self, similarity_data: Dict[str, Any]) -> Any:
        """
        Create a similarity graph from Dashing results.

        Args:
            similarity_data: Dictionary with Dashing results

        Returns:
            Graph representation of similarities
        """
        # Import here to avoid circular imports
        from derep_genomes.core.graph import create_graph_from_edges

        # Extract edges from similarity data
        edges = []
        for i in range(len(similarity_data["source"])):
            source = similarity_data["source"][i]
            target = similarity_data["target"][i]
            weight = similarity_data["weight"][i]
            weight_raw = similarity_data["weight_raw"][i]
            source_len = similarity_data["source_len"][i]
            target_len = similarity_data["target_len"][i]

            # Only add edges between different nodes
            if source != target:
                edges.append(
                    {
                        "source": source,
                        "target": target,
                        "weight": weight,
                        "weight_raw": weight_raw,
                        "aln_frac": 1.0,  # No alignment fraction in Dashing
                        "source_len": source_len,
                        "target_len": target_len,
                    }
                )

        # Create graph
        graph = create_graph_from_edges(edges)

        return graph
