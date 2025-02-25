"""
SKANI implementation for genome similarity calculation.
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


class SKANI(SimilarityMeasure):
    """SKANI-based genome similarity measure."""

    def __init__(
        self,
        threads: int,
        threshold: float,
        temp_dir: Path,
        skani_k: int = 15,
        skani_s: int = 10000,
        skani_min_conf: float = 0.2,
        skani_target_ani: float = 0.95,
        skani_min_kmers: int = 100,
    ):
        """
        Initialize SKANI similarity measure.

        Args:
            threads: Number of threads to use
            threshold: Threshold for similarity (0-1)
            temp_dir: Directory for temporary files
            skani_k: K-mer size for SKANI
            skani_s: Sketch size for SKANI
            skani_min_conf: Minimum confidence value for reporting matches
            skani_target_ani: Target ANI value for optimization
            skani_min_kmers: Minimum number of k-mers two genomes must share
        """
        super().__init__(threads, threshold, temp_dir)
        self.skani_k = skani_k
        self.skani_s = skani_s
        self.skani_min_conf = skani_min_conf
        self.skani_target_ani = skani_target_ani
        self.skani_min_kmers = skani_min_kmers

    def run_skani(self, assembly_list_path: str, output_path: str) -> None:
        """
        Run SKANI on a list of assemblies.

        Args:
            assembly_list_path: Path to file with list of assemblies
            output_path: Path to output file

        Raises:
            subprocess.SubprocessError: If SKANI fails
        """
        cmd = [
            "skani",
            "tri",
            "--list",
            assembly_list_path,
            "-t",
            str(self.threads),
            "-k",
            str(self.skani_k),
            "-s",
            str(self.skani_s),
            "--min-confidence",
            str(self.skani_min_conf),
            "--target-ani",
            str(self.skani_target_ani),
            "--min-shared-kmers",
            str(self.skani_min_kmers),
            "-o",
            output_path,
        ]

        try:
            logger.debug(f"Running SKANI with command: {' '.join(cmd)}")
            subprocess.run(
                cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )
        except subprocess.SubprocessError as e:
            logger.error(f"SKANI failed: {e}")
            raise

    def compute_pairwise(self, assemblies: List[str]) -> Dict[str, Any]:
        """
        Compute pairwise SKANI similarities between assemblies.

        Args:
            assemblies: List of assembly file paths

        Returns:
            Dictionary with pairwise similarity data
        """
        # Prefetch all genome lengths at once
        genome_length_cache.prefetch(assemblies)

        # Create temporary files
        assembly_list_path = self.write_assembly_list(assemblies)
        output_path, _ = self.create_temp_file(suffix=".skani")

        try:
            # Run SKANI
            self.run_skani(assembly_list_path, output_path)

            # Create graph directly from SKANI output
            from derep_genomes.core.graph import create_graph_from_skani

            graph = create_graph_from_skani(output_path)

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
            "ANI": [],
            "af1": [],
            "af2": [],
            "aln_frac": [],
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
            results["aln_frac"].append(data["aln_frac"])
            results["ANI"].append(
                data["weight_raw"] * 100.0
            )  # Convert back to 0-100 scale

            # Infer af1 and af2 from aln_frac (approximately)
            results["af1"].append(data["aln_frac"])
            results["af2"].append(data["aln_frac"])

            results["source_len"].append(graph.nodes[u]["genome_len"])
            results["target_len"].append(graph.nodes[v]["genome_len"])

        return results

    def filter_by_similarity(self, similarity_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Filter SKANI results by similarity threshold.

        Args:
            similarity_data: Dictionary with SKANI results

        Returns:
            Filtered similarity data
        """
        # Convert threshold from 0-1 to ANI scale (0-100)
        ani_threshold = self.threshold * 100.0

        # Find indices of results above threshold
        indices = [
            i for i, ani in enumerate(similarity_data["ANI"]) if ani >= ani_threshold
        ]

        # Filter all keys by indices
        filtered_data = {
            key: [similarity_data[key][i] for i in indices] for key in similarity_data
        }

        logger.debug(
            f"Filtered from {len(similarity_data['ANI'])} to {len(filtered_data['ANI'])} "
            f"comparisons using threshold {self.threshold}"
        )

        return filtered_data

    def create_similarity_graph(self, similarity_data: Dict[str, Any]) -> Any:
        """
        Create a similarity graph from SKANI results.

        Args:
            similarity_data: Dictionary with SKANI results

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
            aln_frac = similarity_data["aln_frac"][i]
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
                        "aln_frac": aln_frac,
                        "source_len": source_len,
                        "target_len": target_len,
                    }
                )

        # Create graph
        graph = create_graph_from_edges(edges)

        return graph
