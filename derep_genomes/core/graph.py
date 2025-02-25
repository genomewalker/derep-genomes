"""
Graph-based clustering algorithms for genome dereplication.
"""

import math
import logging
import random
import os
from typing import List, Dict, Tuple, Set, Union, Optional, Any, Iterator

import networkx as nx
import igraph as ig
import leidenalg
import numpy as np
from scipy import stats
import duckdb

from derep_genomes.utils.helpers import genome_length_cache, get_assembly_length, sanitize_for_filename

logger = logging.getLogger(__name__)


def create_graph_from_edges(edges: List[Dict[str, Any]]) -> nx.Graph:
    """
    Create a NetworkX graph from a list of edges.
    
    Args:
        edges: List of edge dictionaries
        
    Returns:
        NetworkX graph
    """
    graph = nx.Graph()
    
    # Collect all unique nodes
    all_nodes = set()
    for edge in edges:
        all_nodes.add(edge['source'])
        all_nodes.add(edge['target'])
    
    # Prefetch lengths for all nodes at once
    genome_length_cache.prefetch(list(all_nodes))
    
    # Add nodes with genome lengths
    for node in all_nodes:
        graph.add_node(node, genome_len=genome_length_cache.get(node))
    
    # Add edges
    for edge in edges:
        source = edge['source']
        target = edge['target']
        
        # Skip self-loops
        if source == target:
            continue
            
        attrs = {k: v for k, v in edge.items() if k not in ('source', 'target')}
        graph.add_edge(source, target, **attrs)
    
    return graph


def create_graph_from_file(file_path: str, file_type: str = None, threshold: float = 0.01) -> nx.Graph:
    """
    Create graph from similarity/distance file using the most efficient method.
    
    Args:
        file_path: Path to similarity/distance file
        file_type: Type of file ('fastani', 'skani', 'mash', 'dashing', or None for auto-detect)
        threshold: Threshold for filtering (meaning depends on file type)
        
    Returns:
        NetworkX graph
    """
    # Auto-detect file type if not specified
    if file_type is None:
        with open(file_path, 'r') as f:
            first_line = f.readline()
            if '\t' in first_line:
                parts = first_line.strip().split('\t')
                # Check number of columns
                if len(parts) == 5:
                    # Likely FastANI (source, target, ANI, frags, total_frags)
                    file_type = 'fastani'
                elif 'query' in first_line or 'reference' in first_line:
                    # Likely SKANI with header
                    file_type = 'skani'
                else:
                    # Likely distance matrix with header
                    file_type = 'matrix'
    
    # Use the appropriate parser
    if file_type == 'fastani':
        return create_graph_from_fastani(file_path)
    elif file_type == 'skani':
        return create_graph_from_skani(file_path)
    elif file_type in ('mash', 'dashing', 'matrix'):
        return create_graph_from_distance_matrix(file_path, threshold)
    else:
        # Fall back to a generic parser
        return create_graph_generic(file_path)


def create_graph_from_fastani(file_path: str) -> nx.Graph:
    """
    Create graph directly from FastANI output using DuckDB.
    
    Args:
        file_path: Path to FastANI output file
        
    Returns:
        NetworkX graph
    """
    # Create a temporary in-memory DuckDB database
    con = duckdb.connect(":memory:")
    
    # Read the FastANI file format (query, reference, ANI, fragments, total_fragments)
    con.execute(f"""
        CREATE TABLE fastani AS 
        SELECT 
            column0 as source,
            column1 as target,
            CAST(column2 AS FLOAT) as ANI,
            CAST(column3 AS INTEGER) as frags,
            CAST(column4 AS INTEGER) as total_frags
        FROM read_csv('{file_path}', sep='\t', header=False)
    """)
    
    # Calculate additional metrics and extract data
    con.execute("""
        SELECT 
            source, 
            target, 
            ANI / 100.0 AS weight_raw,
            frags / total_frags AS aln_frac,
            (ANI / 100.0) * (frags / total_frags) AS weight
        FROM fastani
        WHERE source != target
    """)
    
    results = con.fetchall()
    
    # Collect all unique genome paths
    all_genomes = set()
    for row in results:
        source, target, _, _, _ = row
        all_genomes.add(source)
        all_genomes.add(target)
    
    # Prefetch genome lengths
    genome_length_cache.prefetch(list(all_genomes))
    
    # Build graph
    graph = nx.Graph()
    
    # Add nodes with genome lengths
    for genome in all_genomes:
        graph.add_node(genome, genome_len=genome_length_cache.get(genome))
    
    # Add edges
    for row in results:
        source, target, weight_raw, aln_frac, weight = row
        graph.add_edge(source, target, 
                      weight=weight, 
                      weight_raw=weight_raw, 
                      aln_frac=aln_frac)
    
    return graph


def create_graph_from_skani(file_path: str) -> nx.Graph:
    """
    Create graph directly from SKANI output using DuckDB.
    
    Args:
        file_path: Path to SKANI output file
        
    Returns:
        NetworkX graph
    """
    # Create a temporary in-memory DuckDB database
    con = duckdb.connect(":memory:")
    
    # Read the SKANI file format
    try:
        # Try to read with header first
        con.execute(f"CREATE TABLE skani AS SELECT * FROM read_csv('{file_path}', sep='\t', header=True)")
        
        # Get the column names
        columns = con.execute("PRAGMA table_info(skani)").fetchall()
        column_names = [col[1] for col in columns]
        
        # SKANI format typically has columns like 'query', 'reference', 'ani', 'af1', 'af2'
        # Try to map column names to expected fields
        source_col = next((col for col in column_names if col.lower() in ('query', 'source')), column_names[0])
        target_col = next((col for col in column_names if col.lower() in ('reference', 'target')), column_names[1])
        ani_col = next((col for col in column_names if col.lower() in ('ani', 'identity')), column_names[2])
        af1_col = next((col for col in column_names if col.lower() in ('af1', 'coverage1')), column_names[3])
        af2_col = next((col for col in column_names if col.lower() in ('af2', 'coverage2')), column_names[4])
        
        # Calculate metrics and extract data
        query = f"""
            SELECT 
                "{source_col}" AS source, 
                "{target_col}" AS target, 
                "{ani_col}" AS weight_raw,
                ("{af1_col}" + "{af2_col}") / 2.0 AS aln_frac,
                "{ani_col}" * (("{af1_col}" + "{af2_col}") / 2.0) AS weight
            FROM skani
            WHERE "{source_col}" != "{target_col}"
        """
        
    except Exception as e:
        # If header parsing fails, try without header
        logger.warning(f"Error parsing SKANI header, trying without header: {e}")
        con.execute("DROP TABLE IF EXISTS skani")
        con.execute(f"""
            CREATE TABLE skani AS 
            SELECT 
                column0 as source,
                column1 as target,
                CAST(column2 AS FLOAT) as ani,
                CAST(column3 AS FLOAT) as af1,
                CAST(column4 AS FLOAT) as af2
            FROM read_csv('{file_path}', sep='\t', header=False)
        """)
        
        query = """
            SELECT 
                source, 
                target, 
                ani AS weight_raw,
                (af1 + af2) / 2.0 AS aln_frac,
                ani * ((af1 + af2) / 2.0) AS weight
            FROM skani
            WHERE source != target
        """
    
    # Execute the query
    con.execute(query)
    results = con.fetchall()
    
    # Collect all unique genome paths
    all_genomes = set()
    for row in results:
        source, target, _, _, _ = row
        all_genomes.add(source)
        all_genomes.add(target)
    
    # Prefetch genome lengths
    genome_length_cache.prefetch(list(all_genomes))
    
    # Build graph
    graph = nx.Graph()
    
    # Add nodes with genome lengths
    for genome in all_genomes:
        graph.add_node(genome, genome_len=genome_length_cache.get(genome))
    
    # Add edges
    for row in results:
        source, target, weight_raw, aln_frac, weight = row
        graph.add_edge(source, target, 
                      weight=weight, 
                      weight_raw=weight_raw, 
                      aln_frac=aln_frac)
    
    return graph


def create_graph_from_distance_matrix(file_path: str, threshold: float = 0.01) -> nx.Graph:
    """
    Create graph from Mash/Dashing distance matrix using DuckDB.
    
    Args:
        file_path: Path to distance matrix file
        threshold: Maximum distance to consider (default: 0.01)
        
    Returns:
        NetworkX graph
    """
    # Create a temporary in-memory DuckDB database
    con = duckdb.connect(":memory:")
    
    # Try to import the matrix as a CSV
    try:
        # First we need to get the header to know the genome names
        with open(file_path, 'r') as f:
            header = f.readline().strip().split('\t')
            genome_names = [h.strip() for h in header[1:]]  # First column is empty/row names
        
        # Now process file with DuckDB
        # Create a temporary wide table for the matrix
        con.execute(f"CREATE TABLE dist_matrix AS SELECT * FROM read_csv('{file_path}', sep='\t', header=True)")
        
        # Determine row names column
        con.execute("PRAGMA table_info(dist_matrix)")
        columns = con.fetchall()
        row_col = columns[0][1]
        
        # Need to transform wide format to long format
        # This is tricky with DuckDB, we'll use a loop to build the query
        unpivot_clauses = []
        for i, genome in enumerate(genome_names):
            col_name = columns[i+1][1]  # +1 to skip row name column
            unpivot_clauses.append(f'"{col_name}" as "{genome}"')
        
        unpivot_query = f"""
            SELECT 
                "{row_col}" as source,
                UNNEST([{', '.join([f"'{g}'" for g in genome_names])}]) as target,
                UNNEST([{', '.join(unpivot_clauses)}]) as distance
            FROM dist_matrix
        """
        
        con.execute(unpivot_query)
        
    except Exception as e:
        logger.error(f"Error processing distance matrix as CSV: {e}")
        logger.info("Falling back to streaming parser for distance matrix")
        
        # Fall back to streaming approach
        return create_graph_from_distance_matrix_streaming(file_path, threshold)
    
    # Filter by threshold and convert to similarity
    con.execute(f"""
        SELECT 
            source, 
            target, 
            CAST(distance AS FLOAT) as distance,
            1.0 - CAST(distance AS FLOAT) as weight,
            1.0 - CAST(distance AS FLOAT) as weight_raw
        FROM ({unpivot_query})
        WHERE distance <= {threshold} AND source != target
    """)
    
    results = con.fetchall()
    
    # Collect all unique genome paths
    all_genomes = set()
    for row in results:
        source, target, _, _, _ = row
        all_genomes.add(source)
        all_genomes.add(target)
    
    # Prefetch genome lengths
    genome_length_cache.prefetch(list(all_genomes))
    
    # Build graph
    graph = nx.Graph()
    
    # Add nodes with genome lengths
    for genome in all_genomes:
        graph.add_node(genome, genome_len=genome_length_cache.get(genome))
    
    # Add edges
    for row in results:
        source, target, distance, weight, weight_raw = row
        graph.add_edge(source, target, 
                      distance=distance,
                      weight=weight, 
                      weight_raw=weight_raw, 
                      aln_frac=1.0)  # No alignment fraction in Mash/Dashing
    
    return graph


def create_graph_from_distance_matrix_streaming(file_path: str, threshold: float = 0.01) -> nx.Graph:
    """
    Create graph from distance matrix using streaming approach for very large files.
    
    Args:
        file_path: Path to distance matrix file
        threshold: Maximum distance to consider
        
    Returns:
        NetworkX graph
    """
    graph = nx.Graph()
    
    # First pass to get all genome names (column/row headers)
    genome_names = []
    with open(file_path, 'r') as f:
        header = f.readline().strip().split('\t')
        genome_names = [h.strip() for h in header[1:]]  # First column is empty/row names
    
    # Prefetch genome lengths
    genome_length_cache.prefetch(genome_names)
    
    # Add nodes
    for genome in genome_names:
        graph.add_node(genome, genome_len=genome_length_cache.get(genome))
    
    # Second pass to read the matrix in chunks
    chunk_size = 1000  # Process 1000 lines at a time
    with open(file_path, 'r') as f:
        # Skip header
        next(f)
        
        for line in f:
            fields = line.strip().split('\t')
            source = fields[0].strip()
            
            # Process distances for this source
            for j, val_str in enumerate(fields[1:]):
                if val_str and val_str != "NA":
                    target = genome_names[j]
                    
                    if source != target:  # Avoid self-loops
                        try:
                            distance = float(val_str)
                            
                            # Only add edges for distances below threshold
                            if distance <= threshold:
                                # Convert distance to similarity (1-distance)
                                weight = 1.0 - distance
                                weight_raw = weight
                                
                                # Add to graph
                                graph.add_edge(source, target, 
                                              distance=distance,
                                              weight=weight, 
                                              weight_raw=weight_raw,
                                              aln_frac=1.0)  # No alignment fraction in Mash/Dashing
                        except ValueError:
                            continue  # Skip non-numeric values
    
    return graph


def create_graph_generic(file_path: str) -> nx.Graph:
    """
    Generic graph creation that tries to handle any format.
    
    Args:
        file_path: Path to similarity/distance file
        
    Returns:
        NetworkX graph
    """
    # Try different parsers in order of likelihood
    try:
        return create_graph_from_fastani(file_path)
    except Exception:
        try:
            return create_graph_from_skani(file_path)
        except Exception:
            try:
                return create_graph_from_distance_matrix(file_path)
            except Exception as e:
                logger.error(f"Failed to parse file {file_path} with any parser: {e}")
                raise ValueError(f"Could not parse file {file_path} with any available parser")


def binary_search_filter(graph: nx.Graph) -> nx.Graph:
    """
    Filter graph by finding the smallest weight where the graph gets disconnected.
    
    Args:
        graph: NetworkX graph
        
    Returns:
        Filtered NetworkX graph
    """
    # Get sorted list of unique edge weights
    weights = sorted(set(data['weight'] for _, _, data in graph.edges(data=True)))
    
    if len(weights) <= 1:
        return graph
    
    g_copy = graph.copy()
    low = 0
    high = len(weights) - 1
    
    while low <= high:
        mid = (low + high) // 2
        threshold = weights[mid]
        
        # Remove edges below threshold
        edges_to_remove = [
            (u, v) for u, v, d in graph.edges(data=True) 
            if d['weight'] <= threshold
        ]
        
        g_test = graph.copy()
        g_test.remove_edges_from(edges_to_remove)
        
        logger.debug(
            f"Filtering -> low: {low} mid: {mid} high: {high} "
            f"threshold: {threshold:.4f} components: {nx.number_connected_components(g_test)} "
            f"edges: {g_test.number_of_edges()}"
        )
        
        if nx.is_connected(g_test):
            # We can try a higher threshold (remove more edges)
            low = mid + 1
            g_copy = g_test
        else:
            # We need a lower threshold (keep more edges)
            high = mid - 1
    
    return g_copy


class Representative:
    """Class for managing representative genomes and their relations."""

    def __init__(self, name: str, graph: nx.Graph, threshold: float = 2.0):
        """
        Initialize the representative genome.
        
        Args:
            name: Name of the representative
            graph: NetworkX graph
            threshold: Z-score threshold for filtering
        """
        self.name = name
        self.n_nodes = len(list(graph.nodes))
        self.graph_avg_weight = None
        self.graph_sd_weight = None
        self.graph_avg_weight_raw = None
        self.graph_sd_weight_raw = None

        self.subgraph_selected_avg_weight = None
        self.subgraph_selected_sd_weight = None
        self.subgraph_selected_avg_weight_raw = None
        self.subgraph_selected_sd_weight_raw = None

        self.subgraph_discarded_avg_weight = None
        self.subgraph_discarded_sd_weight = None
        self.subgraph_discarded_avg_weight_raw = None
        self.subgraph_discarded_sd_weight_raw = None

        # Find genomes similar to this representative
        self.selected = self._refine_candidates(subgraph=graph, threshold=threshold)

        # Calculate global stats for the graph
        (
            self.graph_avg_weight,
            self.graph_sd_weight,
            self.graph_avg_weight_raw,
            self.graph_sd_weight_raw,
        ) = self._get_global_stats(graph)

        # Calculate stats for selected subgraph
        selected = self.selected + [self.name]
        if len(selected) > 1:
            (
                self.subgraph_selected_avg_weight,
                self.subgraph_selected_sd_weight,
                self.subgraph_selected_avg_weight_raw,
                self.subgraph_selected_sd_weight_raw,
            ) = self._get_subgraph_selected(graph, selected)

        # Calculate stats for discarded subgraph
        self.discarded = [x for x in list(graph.nodes) if x not in self.selected and x != self.name]
        if len(self.discarded) > 1:
            (
                self.subgraph_discarded_avg_weight,
                self.subgraph_discarded_sd_weight,
                self.subgraph_discarded_avg_weight_raw,
                self.subgraph_discarded_sd_weight_raw,
            ) = self._get_subgraph_discarded(graph, self.discarded)

        self.n_nodes_selected = len(self.selected)
        self.n_nodes_discarded = len(self.discarded)

    def _get_weight_stats(self, g: nx.Graph) -> Tuple[float, float, float, float]:
        """
        Get weight statistics for a graph.

        Args:
            g: NetworkX graph

        Returns:
            Tuple of (mean_weight, sd_weight, mean_weight_raw, sd_weight_raw)
        """
        weights = []
        weights_raw = []

        for u, v, d in g.edges(data=True):
            weights.append(d["weight"])
            weights_raw.append(d["weight_raw"])

        if len(list(g.edges)) > 1:
            mean_weight = np.mean(weights)
            mean_weight_raw = np.mean(weights_raw)
            sd_weight = np.std(weights)
            sd_weight_raw = np.std(weights_raw)
        else:
            mean_weight = None
            mean_weight_raw = None
            sd_weight = None
            sd_weight_raw = None

        return mean_weight, sd_weight, mean_weight_raw, sd_weight_raw

    def _get_global_stats(self, g: nx.Graph) -> Tuple[float, float, float, float]:
        """
        Get global statistics for a graph.

        Args:
            g: NetworkX graph

        Returns:
            Tuple of (mean_weight, sd_weight, mean_weight_raw, sd_weight_raw)
        """
        return self._get_weight_stats(g)

    def _get_subgraph_selected(
        self, g: nx.Graph, nodes: List[str]
    ) -> Tuple[float, float, float, float]:
        """
        Get statistics for the selected subgraph.

        Args:
            g: NetworkX graph
            nodes: List of node names

        Returns:
            Tuple of (mean_weight, sd_weight, mean_weight_raw, sd_weight_raw)
        """
        subgraph = g.subgraph(nodes)
        return self._get_weight_stats(subgraph)

    def _get_subgraph_discarded(
        self, g: nx.Graph, nodes: List[str]
    ) -> Tuple[float, float, float, float]:
        """
        Get statistics for the discarded subgraph.

        Args:
            g: NetworkX graph
            nodes: List of node names

        Returns:
            Tuple of (mean_weight, sd_weight, mean_weight_raw, sd_weight_raw)
        """
        subgraph = g.subgraph(nodes)
        return self._get_weight_stats(subgraph)

    def _refine_candidates(self, subgraph: nx.Graph, threshold: float) -> List[str]:
        """
        Refine candidate genomes by z-score filtering.

        Args:
            subgraph: NetworkX graph
            threshold: Z-score threshold

        Returns:
            List of selected genome names
        """
        rep = self.name
        results = []

        if subgraph.number_of_edges() == 0:
            return []

        # Collect data for all neighbors
        for neighbor in nx.all_neighbors(subgraph, rep):
            if neighbor != rep:
                source_len = subgraph.nodes[rep]["genome_len"]
                target_len = subgraph.nodes[neighbor]["genome_len"]
                len_diff = abs(source_len - target_len)
                weight = subgraph[rep][neighbor]["weight_raw"]

                # Check if aln_frac is present
                result = {
                    "source": rep,
                    "target": neighbor,
                    "weight": float(weight),
                    "len_diff": int(len_diff),
                }

                if "aln_frac" in subgraph[rep][neighbor]:
                    result["aln_frac"] = float(subgraph[rep][neighbor]["aln_frac"])

                results.append(result)

        # Handle cases with too few results
        if len(results) < 2:
            return []

        # Convert to numpy arrays for efficient z-score calculation
        targets = [r["target"] for r in results]
        weights = np.array([r["weight"] for r in results])
        len_diffs = np.array([r["len_diff"] for r in results])

        has_aln_frac = all("aln_frac" in r for r in results)
        if has_aln_frac:
            aln_fracs = np.array([r["aln_frac"] for r in results])

        # Calculate z-scores
        selected_idx = []

        # Check if values are all the same (avoid division by zero)
        if np.std(weights) > 0:
            weight_zscores = stats.zscore(weights)
            selected_idx.extend(np.where(np.abs(weight_zscores) >= threshold)[0])

        if np.std(len_diffs) > 0:
            len_diff_zscores = stats.zscore(len_diffs)
            selected_idx.extend(np.where(np.abs(len_diff_zscores) >= threshold)[0])

        if has_aln_frac and np.std(aln_fracs) > 0:
            aln_frac_zscores = stats.zscore(aln_fracs)
            selected_idx.extend(np.where(np.abs(aln_frac_zscores) >= threshold)[0])

        # Get unique indices
        selected_idx = list(set(selected_idx))

        # Return selected targets
        return [targets[i] for i in selected_idx]


def is_unique(values: np.ndarray) -> bool:
    """
    Check if all values in an array are identical.

    Args:
        values: NumPy array

    Returns:
        True if all values are identical, False otherwise
    """
    return (values[0] == values).all()


def refine_representatives(
    rep: str, graph: nx.Graph, threshold: float = 2.0
) -> List[str]:
    """
    Refine representative genomes by removing redundant ones.

    Args:
        rep: Name of the representative genome
        graph: NetworkX graph
        threshold: Z-score threshold for filtering

    Returns:
        List of refined representative genomes
    """
    representative_obj = Representative(name=rep, graph=graph, threshold=threshold)
    return representative_obj.selected


def from_networkx_to_igraph(graph: nx.Graph) -> ig.Graph:
    """
    Convert a NetworkX graph to an igraph Graph.

    Args:
        graph: NetworkX graph

    Returns:
        igraph Graph
    """
    # Create igraph Graph
    directed = graph.is_directed()
    ig_graph = ig.Graph(directed=directed)

    # Add vertices with names
    node_names = list(graph.nodes())
    ig_graph.add_vertices(node_names)
    ig_graph.vs["name"] = node_names

    # Add node attributes
    for node, data in graph.nodes(data=True):
        vertex_idx = ig_graph.vs.find(name=node).index
        for key, value in data.items():
            ig_graph.vs[vertex_idx][key] = value

    # Add edges
    edges = [(u, v) for u, v in graph.edges()]
    if edges:  # Only add edges if there are any
        ig_graph.add_edges(edges)

        # Add edge attributes
        for u, v, data in graph.edges(data=True):
            try:
                edge_idx = ig_graph.get_eid(u, v)
                for key, value in data.items():
                    ig_graph.es[edge_idx][key] = value
            except:
                # Edge might not exist in igraph
                pass

    return ig_graph


def get_eigen(g: ig.Graph) -> str:
    """
    Get the node with the highest eigenvector centrality.

    Args:
        g: igraph Graph

    Returns:
        Name of the node with highest centrality
    """
    if g.ecount() == 0:
        if g.vcount() == 1:
            logger.debug("No edges found, one node found")
            return g.vs[0]["name"]
        else:
            logger.debug("No edges found, multiple nodes found")
            return random.choice(g.vs["name"])
    else:
        try:
            cents = g.eigenvector_centrality(directed=False, weights="weight")
            cents = {g.vs[i]["name"]: cents[i] for i in range(len(cents))}
            rep = max(cents.keys(), key=(lambda k: cents[k]))
            return rep
        except:
            # Fallback if eigenvector centrality fails
            logger.debug("Eigenvector centrality failed, using degree centrality")
            cents = g.degree()
            cents = {g.vs[i]["name"]: cents[i] for i in range(len(cents))}
            rep = max(cents.keys(), key=(lambda k: cents[k]))
            return rep


def get_communities(graph: nx.Graph) -> Dict[str, int]:
    """
    Find communities in a graph using the Leiden algorithm.

    Args:
        graph: NetworkX graph

    Returns:
        Dictionary mapping nodes to community IDs
    """
    # Convert NetworkX graph to igraph for Leiden algorithm
    ig_graph = from_networkx_to_igraph(graph)

    # Run Leiden community detection
    logger.debug(
        "Finding genome communities using Leiden community detection algorithm"
    )
    partition = leidenalg.find_partition(
        ig_graph,
        leidenalg.ModularityVertexPartition,
        weights="weight",
        seed=1234,
        n_iterations=-1,
    )

    # Convert partition to dictionary
    communities = [ig_graph.vs[x]["name"] for x in partition]
    partition_dict = {k: i for i, sub_l in enumerate(communities) for k in sub_l}

    return partition_dict


def find_clusters(graph: nx.Graph, threshold: float = 2.0) -> Dict[str, Any]:
    """
    Find clusters in a graph using the Leiden algorithm.

    Args:
        graph: NetworkX graph
        threshold: Z-score threshold for filtering

    Returns:
        Dictionary with cluster information
    """
    if graph.number_of_edges() == 0:
        logger.debug("Graph has no edges, each node is its own cluster")
        return {
            "partition": {node: i for i, node in enumerate(graph.nodes())},
            "representatives": list(graph.nodes()),
            "filtered_representatives": list(graph.nodes()),
        }

    # Filter graph by weight if necessary
    filtered_graph = binary_search_filter(graph)

    # Get communities
    partition = get_communities(filtered_graph)

    # Find representative for each community
    representatives = []
    filtered_representatives = []

    for i in set(partition.values()):
        # Get nodes in this community
        community_nodes = [k for k, v in partition.items() if v == i]

        # Create subgraph for this community
        subgraph = filtered_graph.subgraph(community_nodes)

        # Convert to igraph for eigenvector centrality
        ig_subgraph = from_networkx_to_igraph(subgraph)

        # Find most central node (representative)
        rep = get_eigen(ig_subgraph)
        representatives.append(rep)

        # Refine selection to remove redundant representatives
        filtered_reps = refine_representatives(rep, subgraph, threshold)
        filtered_representatives.append(rep)  # Always include the main rep
        filtered_representatives.extend(filtered_reps)

    logger.debug(
        f"Graph properties: edges={filtered_graph.number_of_edges()}, "
        f"density={nx.density(filtered_graph)}, "
        f"components={nx.number_connected_components(filtered_graph)}, "
        f"communities={len(set(partition.values()))}"
    )

    return {
        "partition": partition,
        "representatives": representatives,
        "filtered_representatives": filtered_representatives,
    }
