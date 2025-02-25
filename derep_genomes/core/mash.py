"""
Mash implementation for genome similarity pre-clustering.
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


class Mash(SimilarityMeasure):
    """Mash-based genome similarity measure for pre-clustering."""

    def __init__(
        self,
        threads: int,
        threshold: float,
        temp_dir: Path,
        mash_k: int = 21,
        mash_s: int = 10000,
        mash_threshold: float = 0.01,
    ):
        """
        Initialize Mash similarity measure.

        Args:
            threads: Number of threads to use
            threshold: Threshold for similarity (0-1)
            temp_dir: Directory for temporary files
            mash_k: K-mer size for Mash
            mash_s: Sketch size for Mash
            mash_threshold: Maximum Mash distance to consider genomes as similar
        """
        super().__init__(threads, threshold, temp_dir)
        self.mash_k = mash_k
        self.mash_s = mash_s
        self.mash_threshold = mash_threshold

    def create_mash_sketch(self, assembly_list_path: str) -> str:
        """
        Create Mash sketch from a list of assemblies.

        Args:
            assembly_list_path: Path to file with list of assemblies

        Returns:
            Path to Mash sketch file
        """
        sketch_path = os.path.join(self.temp_dir, "mash.msh")

        cmd = [
            "mash",
            "sketch",
            "-p",
            str(self.threads),
            "-o",
            os.path.join(self.temp_dir, "mash"),
            "-s",
            str(self.mash_s),
            "-k",
            str(self.mash_k),
            "-l",
            assembly_list_path,
        ]

        try:
            logger.debug(f"Creating Mash sketch with command: {' '.join(cmd)}")
            subprocess.run(
                cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )
        except subprocess.SubprocessError as e:
            logger.error(f"Mash sketch creation failed: {e}")
            raise

        return sketch_path

    def run_mash(self, sketch_path: str, output_path: str) -> None:
        """
        Run Mash on a sketch file.

        Args:
            sketch_path: Path to Mash sketch file
            output_path: Path to output file

        Raises:
            subprocess.SubprocessError: If Mash fails
        """
        cmd = ["mash", "dist", "-p", str(self.threads), sketch_path, sketch_path]

        try:
            logger.debug(f"Running Mash with command: {' '.join(cmd)}")
            with open(output_path, "w") as f:
                subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE)
        except subprocess.SubprocessError as e:
            logger.error(f"Mash failed: {e}")
            raise

    def compute_pairwise(self, assemblies: List[str]) -> Dict[str, Any]:
        """
        Compute pairwise Mash similarities between assemblies.

        Args:
            assemblies: List of assembly file paths

        Returns:
            Dictionary with pairwise similarity data
        """
        # Prefetch all genome lengths at once
        genome_length_cache.prefetch(assemblies)

        # Create temporary files
        assembly_list_path = self.write_assembly_list(assemblies)
        output_path, _ = self.create_temp_file(suffix=".mash")

        try:
            # Create Mash sketch
            sketch_path = self.create_mash_sketch(assembly_list_path)

            # Run Mash
            self.run_mash(sketch_path, output_path)

            # Create graph directly from Mash output
            from derep_genomes.core.graph import create_graph_from_distance_matrix

            graph = create_graph_from_distance_matrix(output_path, self.mash_threshold)

            # Convert graph to edge list for consistent interface
            results = self._graph_to_edge_list(graph)

            return results
        finally:
            # Clean up temporary files
            try:
                os.remove(assembly_list_path)
                os.remove(os.path.join(self.temp_dir, "mash.msh"))
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
        Filter Mash results by distance threshold.

        Args:
            similarity_data: Dictionary with Mash results

        Returns:
            Filtered similarity data
        """
        # Find indices of results where distance is below threshold
        indices = [
            i
            for i, distance in enumerate(similarity_data.get("distance", []))
            if distance <= self.mash_threshold
        ]

        # Filter all keys by indices
        filtered_data = {
            key: [similarity_data[key][i] for i in indices] for key in similarity_data
        }

        logger.debug(
            f"Filtered from {len(similarity_data.get('distance', []))} to {len(filtered_data.get('distance', []))} "
            f"comparisons using distance threshold {self.mash_threshold}"
        )

        return filtered_data

    def create_similarity_graph(self, similarity_data: Dict[str, Any]) -> Any:
        """
        Create a similarity graph from Mash results.

        Args:
            similarity_data: Dictionary with Mash results

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
                        "aln_frac": 1.0,  # No alignment fraction in Mash
                        "source_len": source_len,
                        "target_len": target_len,
                    }
                )

        # Create graph
        graph = create_graph_from_edges(edges)

        return graph
