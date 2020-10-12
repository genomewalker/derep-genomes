import networkx as nx
import math
import tempfile
from statistics import mean
from community import community_louvain
import pandas as pd
from scipy import stats
import subprocess
import pandas as pd
from derep_genomes.general import get_open_func, get_assembly_length, get_contig_lengths
import logging
from pathlib import Path
import os
from itertools import product
from simple_slurm import Slurm
import yaml
from multiprocessing import Pool
import tqdm
import io, sys

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


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
    return rfile


def process_fastANI_results(rfile):
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


def binary_search_filter(g, low, high, weights, debug):
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
        if debug:
            logging.info(
                "Filtering -> low: {} mid: {} high: {} mid_weight: {} components: {} edges: {}".format(
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
            if debug:
                logging.info(
                    "Final -> low: {} mid: {} high: {} mid_weight: {} components: {} edges: {}".format(
                        low,
                        mid,
                        high,
                        weights[mid],
                        nx.number_connected_components(x),
                        x.number_of_edges(),
                    )
                )
            if nx.number_connected_components(x) == 1:
                return x, weights[mid]
            else:
                logging.info(
                    "Couldn't find a filtering threshold, returning original graph"
                )
                return g, None

        if (high - low) == 0:
            edges2rm = [
                (u, v)
                for u, v, d in g.edges(data=True)
                if float(d["weight"]) <= float(weights[mid + 1])
            ]
            x.remove_edges_from(edges2rm)
            if debug:
                logging.info(
                    "low: {} mid: {} high: {} mid_weight: {} components: {} edges: {}".format(
                        low,
                        mid,
                        high,
                        weights[mid],
                        nx.number_connected_components(x),
                        x.number_of_edges(),
                    )
                )
            if nx.number_connected_components(x) == 1:
                print(nx.number_connected_components(x))
                return x, weights[mid]
            else:
                logging.info(
                    "Couldn't find a filtering threshold, returning original graph"
                )
                return g, None

        if nx.is_connected(x):
            x, w = binary_search_filter(
                g, low=mid, high=high, weights=weights, debug=debug
            )
        else:
            x, w = binary_search_filter(
                g, low=low, high=mid, weights=weights, debug=debug
            )

        if nx.number_connected_components(x) == 1:
            return x, w
        else:
            logging.info(
                "Couldn't find a filtering threshold, returning original graph"
            )
            return g, None


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


def save_chunks_to_disk(chunks, temp_dir):
    wdir = os.path.join(temp_dir, "chunks")
    Path(wdir).mkdir(parents=True, exist_ok=True)
    file_chunks = {}
    for i, assemblies in enumerate(chunks):
        fname = os.path.join(wdir, "chunk_" + str(i))
        with open(fname, "w") as outfile:
            for assm in assemblies:
                outfile.write("%s\n" % assm)
        file_chunks[i] = fname
    return file_chunks, wdir


def create_slurm_commands(files, wdir):
    file_list = list(product(files.keys(), files.keys()))
    odir = os.path.join(wdir, "fastANI_out")
    Path(odir).mkdir(parents=True, exist_ok=True)
    cmds = []
    ofiles = []
    for file in file_list:
        ofile = os.path.join(odir, "fastANI_out_" + str(file[0]) + "_" + str(file[1]))
        cmd = " ".join(
            ["fastANI", "--ql", files[file[0]], "--rl", files[file[1]], "-o", ofile]
        )
        cmds.append(cmd)
        ofiles.append(ofile)
    return cmds, ofiles, odir


def check_slurm_output(files):
    pass


def split_fixed_size(lst, n):
    splits = [lst[i * n : (i + 1) * n] for i in range((len(lst) + n - 1) // n)]
    return splits


def map_slurm_jobs(cmds, ofiles, slurm_config, odir, max_jobs_array):
    rfile = os.path.join(odir, "fastANI_out_jobs.txt")
    with open(rfile, "w") as outfile:
        for cmd in cmds:
            outfile.write("%s\n" % cmd)

    slurm = Slurm(**yaml.load(slurm_config, Loader=yaml.FullLoader))

    if len(cmds) > max_jobs_array:
        logging.info(
            "The number of jobs ({}) is larger than the allowed array size ({}). Splitting jobs..".format(
                len(cmds), max_jobs_array
            )
        )

    job_ranges = split_fixed_size(range(1, len(cmds) + 1), max_jobs_array)
    job_ids = []
    logging.info(
        "Mapping {} jobs to SLURM in {} batch(es)".format(len(cmds), len(job_ranges))
    )
    for r in job_ranges:
        slurm.set_array(r)
        slurm.add_arguments(wait="")
        bjob = "$(awk -vJ=" + "$SLURM_ARRAY_TASK_ID " + "'NR==J' " + rfile + " )"
        text_trap = io.StringIO()
        sys.stdout = text_trap
        job_id = slurm.sbatch(bjob)
        sys.stdout = sys.__stdout__
        job_ids.append(job_id)

    return job_ids


def reduce_slurm_jobs(ofiles, threads):
    if debug is True:
        dfs = list(map(process_fastANI_results, ofiles))
    else:
        p = Pool(threads)
        dfs = list(
            tqdm.tqdm(
                p.imap_unordered(process_fastANI_results, ofiles), total=len(ofiles)
            )
        )
    dfs = pd.concat(dfs)
    return dfs


def dereplicate(
    all_assemblies,
    out_dir,
    threads,
    threshold,
    chunks,
    slurm_config,
    tmp_dir,
    max_jobs_array,
    con,
):
    """
    This function dereplicates genomes by Taxon by:
    1. Calculating ANI and fraction aligned
    2. Dynamically filters the ANI graph
    3.
    """
    derep_assemblies = []
    n_assemblies = all_assemblies.shape[0]
    with tempfile.TemporaryDirectory(dir=tmp_dir, prefix="gderep-") as temp_dir:

        if n_assemblies <= 10 or slurm_config is None:
            if (n_assemblies * n_assemblies) < threads:
                threads = n_assemblies * n_assemblies
            logging.info(
                "Found {} assemblies, using default fastANI with {} threads".format(
                    n_assemblies, threads
                )
            )
            ani_results = pairwise_fastANI(
                assemblies=all_assemblies["assemblies"].tolist(), threads=threads, temp_dir=temp_dir
            )
            logging.info("Processing ANI results")
            pairwise_distances = process_fastANI_results(ani_results)
        else:
            n = chunks
            chunks = split_fixed_size(all_assemblies["assemblies"].tolist(), n)
            # chunks = [
            #     all_assemblies[i * n : (i + 1) * n]
            #     for i in range((len(all_assemblies) + n - 1) // n)
            # ]
            logging.info(
                "Found {} assemblies, creating {} chunks with {} genomes".format(
                    len(all_assemblies), len(chunks), n
                )
            )
            # save chunks to disk
            files, wdir = save_chunks_to_disk(chunks, temp_dir)
            # Create fastANI commands
            cmds, ofiles, odir = create_slurm_commands(files, wdir)
            # Run slurm array job
            slurm_jobs = map_slurm_jobs(
                cmds, ofiles, slurm_config, odir, max_jobs_array
            )
            logging.info("Reducing SLURM jobs and processing ANI results")
            pairwise_distances = reduce_slurm_jobs(ofiles, threads)
        # Convert pw dist to graph
        logging.info("Generating ANI graph")
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
        logging.info("Filtering ANI graph")
        G_filt, w_filt = binary_search_filter(
            g=G, low=0, high=len(weights) - 1, weights=weights, debug=debug
        )

        logging.info(
            "Finding genome representatives using Louvain + eigenvector centrality"
        )

        partition = community_louvain.best_partition(G_filt, resolution=1.0)

        logging.info(
            "Graph properties: w_filt={} components={} edges={} communities={}".format(
                w_filt,
                nx.number_connected_components(G_filt),
                G_filt.number_of_edges(),
                len(set([partition[k] for k in partition])),
            )
        )
        logging.info(
            "Refining genome selection (length difference and fraction aligned, z-score={})".format(
                threshold
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
            if len(candidates) > 1:
                for candidate in candidates:
                    derep_assemblies.append(candidate)
            else:
                derep_assemblies = candidates
    derep_assemblies = list(set(derep_assemblies))
    logging.info(
        "Keeping {}/{} genomes".format(len(derep_assemblies), len(all_assemblies))
    )
    results = (
        w_filt,
        len(set([partition[k] for k in partition])),
        len(all_assemblies),
        len(derep_assemblies),
    )
    rep_keys = [k for k in acc_to_assemblies if acc_to_assemblies[k] in reps]
    return derep_assemblies, results, rep_keys


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

    if len(df.index) < 2:
        if (source_len / len_diff) >= 0.1:
            assms = df["target"].tolist()
            assms.append(rep)
        else:
            assms = [rep]
    else:
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
