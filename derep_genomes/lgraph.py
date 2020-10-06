import networkx as nx
import math
import tempfile
from statistics import mean
import community as community_louvain
import pandas as pd
from scipy import stats
import subprocess
import pandas as pd


def pairwise_fastANI(assemblies, threads, temp_dir):
    glist = temp_dir + "/assm_list.txt"
    rfile = temp_dir + "/fastANI.txt"
    with open(glist, "w") as outfile:
        for assm in assemblies:
            outfile.write("%s\n" % assm)

    fastANI_cmd = [
        "fastANI",
        "--ql",
        glist,
        "--rl",
        glist,
        "-t",
        str(threads),
        "-o",
        rfile,
    ]
    subprocess.run(fastANI_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    df = pd.read_csv(
        rfile,
        sep="\t",
        names=["source", "target", "ANI", "frags", "total_frags"],
    )

    df["aln_frac"] = df["frags"] / df["total_frags"]
    df["weight"] = df["ANI"].div(100)

    df1 = pd.DataFrame(
        set(df["source"].tolist() + df["target"].tolist()),
        columns=["assm"],
    )

    df1["len"] = df1["assm"].map(lambda x: get_assembly_length(x))

    df = df.merge(df1.rename(columns={"assm": "source", "len": "source_len"}))
    df = df.merge(df1.rename(columns={"assm": "target", "len": "target_len"}))
    return df


def binary_search_filter(g, low, high, weights):
    """
    This functions filters a graph by finding the smallest weight where
    the graph gets disconnected.
    """
    x = g.copy()
    edges2rm = []

    if high < low:
        return None, None
    else:
        mid = math.floor((low + high) / 2)
        edges2rm = [
            (u, v)
            for u, v, d in g.edges(data=True)
            if float(d["weight"]) <= float(weights[mid])
        ]
        x.remove_edges_from(edges2rm)
        print(
            "low: {} mid: {} high: {} mid_weight: {} components: {} edges: {}".format(
                low,
                mid,
                high,
                weights[mid],
                nx.number_connected_components(x),
                x.number_of_edges(),
            )
        )

        if (mid - low) == 0:
            edges2rm = [
                (u, v)
                for u, v, d in g.edges(data=True)
                if float(d["weight"]) <= float(weights[mid])
            ]
            x.remove_edges_from(edges2rm)
            print(
                "Final -> low: {} mid: {} high: {} mid_weight: {} components: {} edges: {}".format(
                    low,
                    mid,
                    high,
                    weights[mid],
                    nx.number_connected_components(x),
                    x.number_of_edges(),
                )
            )
            return x, weights[mid]

        if (high - low) == 0:
            edges2rm = [
                (u, v)
                for u, v, d in g.edges(data=True)
                if float(d["weight"]) <= float(weights[mid + 1])
            ]
            x.remove_edges_from(edges2rm)
            print(
                "low: {} mid: {} high: {} mid_weight: {} components: {} edges: {}".format(
                    low,
                    mid,
                    high,
                    weights[mid],
                    nx.number_connected_components(x),
                    x.number_of_edges(),
                )
            )
            return x, weights[mid]

        if nx.is_connected(x):
            x = binary_search_filter(g, low=mid, high=high, weights=weights)
        else:
            x = binary_search_filter(g, low=low, high=mid, weights=weights)
        return x


def get_reps(graph, partition):
    """
    This function finds the representative genomes based on centrality measures
    """
    subgraphs = []
    reps = []
    for i in set([partition[k] for k in partition]):
        nodes = [k for k in partition if partition[k] == i]
        subgraph = graph.subgraph(nodes)
        subgraphs.append(subgraph)
        cents = nx.eigenvector_centrality(subgraph, weight="weight")
        rep = max(cents.keys(), key=(lambda k: cents[k]))
        reps.append(rep)
    return reps, subgraphs


def dereplicate(acc_to_assemblies, threads, threshold):
    """
    This function dereplicates genomes by Taxon by:
    1.
    """
    all_assemblies = sorted(acc_to_assemblies.values())
    derep_assemblies = []
    with tempfile.TemporaryDirectory() as temp_dir:
        # mash_sketch = build_mash_sketch(all_assemblies, threads, temp_dir)
        pairwise_distances = pairwise_fastANI(
            assemblies=all_assemblies, threads=threads, temp_dir=temp_dir
        )

        # Convert pw dist to graph
        M = nx.from_pandas_edgelist(
            pairwise_distances, edge_attr=True, create_using=nx.MultiGraph()
        )
        G = nx.Graph()
        for u, v, data in M.edges(data=True):
            if not G.has_edge(u, v):
                # set weight to 1 if no weight is given for edge in M
                weight = mean(
                    d.get("weight", 1) for d in M.get_edge_data(u, v).values()
                )
                G.add_edge(u, v, weight=weight)

        G.remove_edges_from(nx.selfloop_edges(G))

        weights = []
        for u, v, d in G.edges(data=True):
            weights.append(d["weight"])

        weights = sorted(set(weights))
        print("Filtering networks...")
        G_filt, w_filt = binary_search_filter(
            g=G, low=0, high=len(weights) - 1, weights=weights
        )
        partition = community_louvain.best_partition(G_filt, resolution=1.0)

        print(
            "w_filt: {} components: {} edges: {} communities: {}".format(
                w_filt,
                nx.number_connected_components(G_filt),
                G_filt.number_of_edges(),
                len(set([partition[k] for k in partition])),
            )
        )

        reps, subgraphs = get_reps(graph=G_filt, partition=partition)
        derep_assemblies = []
        for i in range(len(reps)):
            candidates = refine_candidates(
                rep=reps[i],
                subgraph=subgraphs[i],
                threshold=threshold,
                pw=pairwise_distances,
            )
            print("ref: {} -> {}".format(reps[i], candidates))
            if len(candidates) > 1:
                for candidate in candidates:
                    derep_assemblies.append(candidate)
            else:
                derep_assemblies.append(candidate)

    return set(derep_assemblies)


def refine_candidates(rep, subgraph, pw, threshold=2.0):
    results = []
    for reachable_node in nx.dfs_postorder_nodes(subgraph, source=rep):
        if reachable_node != rep:
            df = pw[
                ((pw["source"] == rep) & (pw["target"] == reachable_node))
                | ((pw["target"] == rep) & (pw["source"] == reachable_node))
            ]
            source_len = df[df["source"] == rep].iloc[0]["source_len"]
            target_len = df[df["source"] == rep].iloc[0]["target_len"]
            len_diff = abs(source_len - target_len)

            source_aln = df[df["source"] == rep].iloc[0]["aln_frac"]
            target_aln = df[df["source"] != rep].iloc[0]["aln_frac"]
            aln_frac = mean([source_aln, target_aln])

            source_w = df[df["source"] == rep].iloc[0]["weight"]
            target_w = df[df["source"] != rep].iloc[0]["weight"]
            weight = mean([source_w, target_w])

            results.append(
                {
                    "source": rep,
                    "target": reachable_node,
                    "weight": float(weight),
                    "aln_frac": float(aln_frac),
                    "len_diff": int(len_diff),
                }
            )
    df = pd.DataFrame(results)
    df["z_score_aln"] = stats.zscore(df["aln_frac"])
    df["z_score_diff"] = stats.zscore(df["len_diff"])

    df = df[
        (df["z_score_aln"] <= (-1 * threshold))
        | (df["z_score_diff"].abs() >= threshold)
    ]
    assms = df["target"].tolist()
    assms.append(rep)
    if df.empty:
        return [rep]
    else:
        return assms


def get_assembly_length(filename):
    contig_lengths = sorted(get_contig_lengths(filename), reverse=True)
    total_length = sum(contig_lengths)
    return total_length


def get_contig_lengths(filename):
    lengths = []
    with get_open_func(filename)(filename, "rt") as fasta_file:
        name = ""
        sequence = ""
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":  # Header line = start of new contig
                if name:
                    lengths.append(len(sequence))
                    sequence = ""
                name = line[1:].split()[0]
            else:
                sequence += line
        if name:
            lengths.append(len(sequence))
    return lengths