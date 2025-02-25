"""
FastANI implementation for genome similarity calculation.
"""

import os
import math
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


class FastANI(SimilarityMeasure):
    """FastANI-based genome similarity measure."""

    def __init__(
        self,
        threads: int,
        threshold: float,
        temp_dir: Path,
        ani_fraglen: Optional[int] = None,
        ani_fraglen_fraction: float = 0.005,
        ani_min_frac: float = 0.0,
        ani_threshold: float = 95.0,
        min_genome_size: int = 20000000,
    ):
        """
        Initialize FastANI similarity measure.

        Args:
            threads: Number of threads to use
            threshold: Threshold for similarity (0-1)
            temp_dir: Directory for temporary files
            ani_fraglen: Fragment length for FastANI (None for auto-calculation)
            ani_fraglen_fraction: Fraction of genome size for fragment length calculation
            ani_min_frac: Minimum fraction of genome that must align
            ani_threshold: Minimum ANI value to consider genomes as similar
            min_genome_size: Minimum genome size for fragment length calculation
        """
        super().__init__(threads, threshold, temp_dir)
        self.ani_fraglen = ani_fraglen
        self.ani_fraglen_fraction = ani_fraglen_fraction
        self.ani_min_frac = ani_min_frac
        self.ani_threshold = ani_threshold
        self.min_genome_size = min_genome_size

    def estimate_fragment_length(self, assemblies: List[str]) -> int:
        """
        Estimate appropriate fragment length for FastANI.

        Args:
            assemblies: List of assembly file paths

        Returns:
            Estimated fragment length
        """
        # If fragment length is explicitly set, use it
        if self.ani_fraglen is not None:
            logger.debug(f"Using user-specified fragment length: {self.ani_fraglen}")
            return self.ani_fraglen

        # Prefetch all assembly lengths
        genome_length_cache.prefetch(assemblies)

        # Calculate based on genome size
        assembly_lengths = [
            genome_length_cache.get(assembly) for assembly in assemblies
        ]
        min_length = min(assembly_lengths)

        if min_length < 3000:
            # Round down to nearest 500 for very small genomes
            frag_len = math.floor(min_length / 500.0) * 500.0
        else:
            if min_length > self.min_genome_size:
                frag_len = int(self.ani_fraglen_fraction * min_length)
            else:
                frag_len = 3000

        logger.debug(
            f"Minimum assembly length: {min_length}. "
            f"FastANI fragment length used: {frag_len}"
        )

        if frag_len == 0:
            raise ValueError(
                "Estimated fragment length is 0, assemblies may be too short"
            )

        return frag_len

    def run_fastani(
        self, assembly_list_path: str, frag_len: int, output_path: str
    ) -> None:
        """
        Run FastANI on a list of assemblies.

        Args:
            assembly_list_path: Path to file with list of assemblies
            frag_len: Fragment length to use
            output_path: Path to output file

        Raises:
            subprocess.SubprocessError: If FastANI fails
        """
        cmd = [
            "fastANI",
            "--ql",
            assembly_list_path,
            "--rl",
            assembly_list_path,
            "--fragLen",
            str(int(frag_len)),
            "-t",
            str(self.threads),
            "-o",
            output_path,
        ]

        # Add minimum fraction parameter if set
        if self.ani_min_frac > 0:
            cmd.extend(["--minFraction", str(self.ani_min_frac)])

        try:
            logger.debug(f"Running FastANI with command: {' '.join(cmd)}")
            subprocess.run(
                cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )
        except subprocess.SubprocessError as e:
            logger.error(f"FastANI failed: {e}")
            raise

    def compute_pairwise(self, assemblies: List[str]) -> Dict[str, Any]:
        """
        Compute pairwise FastANI similarities between assemblies.

        Args:
            assemblies: List of assembly file paths

        Returns:
            Dictionary with pairwise similarity data
        """
        # Prefetch genome lengths
        genome_length_cache.prefetch(assemblies)

        # Create temporary files
        assembly_list_path = self.write_assembly_list(assemblies)
        output_path, _ = self.create_temp_file(suffix=".fastani")

        try:
            # Estimate fragment length
            frag_len = self.estimate_fragment_length(assemblies)

            # Run FastANI
            self.run_fastani(assembly_list_path, frag_len, output_path)

            # Create graph directly from FastANI output
            from derep_genomes.core.graph import create_graph_from_fastani

            graph = create_graph_from_fastani(output_path)

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
            "frags": [],
            "total_frags": [],
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

            # Infer frags and total_frags from aln_frac
            frags = int(data["aln_frac"] * 1000)  # Approximate
            total_frags = 1000
            results["frags"].append(frags)
            results["total_frags"].append(total_frags)

            results["source_len"].append(graph.nodes[u]["genome_len"])
            results["target_len"].append(graph.nodes[v]["genome_len"])

        return results

    def filter_by_similarity(self, similarity_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Filter FastANI results by similarity threshold.

        Args:
            similarity_data: Dictionary with FastANI results

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
        Create a similarity graph from FastANI results.

        Args:
            similarity_data: Dictionary with FastANI results

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
