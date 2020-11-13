import math
import tempfile
from statistics import mean, median, stdev
import leidenalg
import networkx as nx
import igraph as ig
from networkx.generators.geometric import thresholded_random_geometric_graph
from numpy.lib.function_base import select
import pandas as pd
from scipy import stats
import subprocess
import pandas as pd
from derep_genomes.general import (
    get_open_func,
    get_assembly_length,
    get_contig_lengths,
    get_assembly_n50,
    is_debug,
)
import logging
from pathlib import Path
import os
from itertools import product, chain
from simple_slurm import Slurm
import yaml
from multiprocessing import Pool
import tqdm
import io, sys
from multiprocessing.pool import ThreadPool

log = logging.getLogger("my_logger")


class representative:
    def __init__(self, name, graph, pw, threshold=2.0):
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

        self.selected = self.refine_candidates(
            subgraph=graph, pw=pw, threshold=threshold
        )

        (
            self.graph_avg_weight,
            self.graph_sd_weight,
            self.graph_avg_weight_raw,
            self.graph_sd_weight_raw,
        ) = self.get_global_stats(g=graph)

        selected = self.selected + [self.name]

        if len(selected) > 1:
            (
                self.subgraph_selected_avg_weight,
                self.subgraph_selected_sd_weight,
                self.subgraph_selected_avg_weight_raw,
                self.subgraph_selected_sd_weight_raw,
            ) = self.get_subgraph_selected(g=graph, nodes=selected)

        self.discarded = [x for x in list(graph.nodes) if x not in self.selected]

        discarded = self.discarded

        if len(discarded) > 1:
            (
                self.subgraph_discarded_avg_weight,
                self.subgraph_discarded_sd_weight,
                self.subgraph_discarded_avg_weight_raw,
                self.subgraph_discarded_sd_weight_raw,
            ) = self.get_subgraph_discarded(g=graph, nodes=discarded)

        self.n_nodes_selected = len(self.selected)
        self.n_nodes_discarded = len(self.discarded) - 1

    def get_weight_stats(self, g):
        weights = []
        weights_raw = []
        for u, v, d in g.edges(data=True):
            weights.append(d["weight"])
            weights_raw.append(d["weight_raw"])
        # if len(list(g.edges)) > 1:
        mean_weight = mean(weights)
        mean_weight_raw = mean(weights_raw)
        if len(list(g.edges)) > 1:
            sd_weight = stdev(weights)
            sd_weight_raw = stdev(weights_raw)
        else:
            sd_weight = None
            sd_weight_raw = None

        return mean_weight, sd_weight, mean_weight_raw, sd_weight_raw

    def get_global_stats(self, g):
        mean_weight, sd_weight, mean_weight_raw, sd_weight_raw = self.get_weight_stats(
            g=g
        )
        return mean_weight, sd_weight, mean_weight_raw, sd_weight_raw

    def get_subgraph_selected(self, g, nodes):
        weights = []
        weights_raw = []
        subgraph = nx.subgraph(g, nodes)
        mean_weight, sd_weight, mean_weight_raw, sd_weight_raw = self.get_weight_stats(
            g=subgraph
        )
        return mean_weight, sd_weight, mean_weight_raw, sd_weight_raw

    def get_subgraph_discarded(self, g, nodes):
        weights = []
        weights_raw = []
        subgraph = nx.subgraph(g, nodes)
        mean_weight, sd_weight, mean_weight_raw, sd_weight_raw = self.get_weight_stats(
            g=subgraph
        )
        return mean_weight, sd_weight, mean_weight_raw, sd_weight_raw

    def refine_candidates(self, subgraph, pw, threshold):
        rep = self.name
        results = []
        len_diff = 0
        source_len = 0
        if subgraph.number_of_edges() == 0:
            return [rep]
        for reachable_node in nx.all_neighbors(subgraph, rep):
            if reachable_node != rep:
                idx = pw["target"] == rep
                pw.loc[idx, ["source", "target", "source_len", "target_len"]] = pw.loc[
                    idx, ["target", "source", "target_len", "source_len"]
                ].values

                df = pw[(pw["source"] == rep) & (pw["target"] == reachable_node)]

                # Switch column order
                # idx = df["target"] == rep
                # df.loc[idx, ["source", "target", "source_len", "target_len"]] = df.loc[
                #     idx, ["target", "source", "target_len", "source_len"]
                # ].values

                source_len = mean(df.source_len.values)
                target_len = mean(df.target_len.values)
                len_diff = abs(source_len - target_len)
                weight = mean(df.weight.values)

                if "aln_frac" in df.columns:
                    aln_frac = mean(df.aln_frac.values)
                    results.append(
                        {
                            "source": rep,
                            "target": reachable_node,
                            "weight": float(weight),
                            "aln_frac": float(aln_frac),
                            "len_diff": int(len_diff),
                        }
                    )
                else:
                    results.append(
                        {
                            "source": rep,
                            "target": reachable_node,
                            "weight": float(weight),
                            "len_diff": int(len_diff),
                        }
                    )
        df = pd.DataFrame(results)

        if len(df.index) < 2:
            if len_diff > 0:
                diff_ratio = source_len / len_diff
            else:
                diff_ratio = 0

            if diff_ratio >= 0.1:
                assms = df["target"].tolist()
                # assms.append(rep)
            else:
                assms = []
        else:
            if "aln_frac" in df.columns:
                if is_unique(df["aln_frac"]):
                    df["z_score_aln"] = 0
                else:
                    df["z_score_aln"] = stats.zscore(df["aln_frac"])

            if is_unique(df["len_diff"]):
                df["z_score_diff"] = 0
            else:
                df["z_score_diff"] = stats.zscore(df["len_diff"])

            if "aln_frac" in df.columns:
                df = df[
                    (df["z_score_aln"].abs() >= threshold)
                    | (df["z_score_diff"].abs() >= threshold)
                ]
            else:
                df = df[(df["z_score_diff"].abs() >= threshold)]

            assms = df["target"].tolist()

        if df.empty:
            return []
        else:
            return assms


# From https://github.com/GiulioRossetti/cdlib/blob/master/cdlib/utils.py#L79
def __from_nx_to_igraph(g, directed=None):
    """
    :param g:
    :param directed:
    :return:
    """

    if ig is None:
        raise ModuleNotFoundError(
            "Optional dependency not satisfied: install igraph to use the selected feature."
        )

    if directed is None:
        directed = g.is_directed()

    gi = ig.Graph(directed=directed)

    ## Two problems to handle:
    # 1)in igraph, names have to be str.
    # 2)since we can ask to compute metrics with found communities and the the original graph, we need to keep
    # the original nodes types in communities. Therefore we need to handle some transparent conversion for non-str nodes
    if type(list(g.nodes)[0]) is str:  # if nodes are string, no problem
        gi.add_vertices([n for n in g.nodes()])
        gi.add_edges([(u, v) for (u, v) in g.edges()])

    else:
        if set(range(len(g.nodes))) == set(
            g.nodes()
        ):  # if original names are well formed contiguous ints, keep this for efficiency.
            # Put these int as str with identitiers in the name attribute
            gi.add_vertices(len(g.nodes))
            gi.add_edges([(u, v) for (u, v) in g.edges()])
            gi.vs["name"] = ["\\" + str(n) for n in g.nodes()]
        else:  # if names are not well formed ints, convert to string and use the identifier to remember
            # converting back to int
            # convert = {str(x):x for x in g.nodes()}
            gi.add_vertices(["\\" + str(n) for n in g.nodes()])
            gi.add_edges([("\\" + str(u), "\\" + str(v)) for (u, v) in g.edges()])

    edgelist = nx.to_pandas_edgelist(g)
    for attr in edgelist.columns[2:]:
        gi.es[attr] = edgelist[attr]

    return gi


def pairwise_fastANI(assemblies, threads, temp_dir, frag_len):
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
        "--fragLen",
        str(int(frag_len)),
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
    # log.debug("Loaded {} comparinsons from ANI file".format(df.shape[0]))
    df["aln_frac"] = df["frags"] / df["total_frags"]
    df["weight"] = df["ANI"].div(100)
    df["weight_raw"] = df["ANI"].div(100)
    df["weight"] = df["weight"] * df["aln_frac"]
    df1 = pd.DataFrame(
        set(df["source"].tolist() + df["target"].tolist()),
        columns=["assm"],
    )

    df1["len"] = df1["assm"].map(lambda x: get_assembly_length(x))
    # df1["n50"] = df1["assm"].map(lambda x: get_assembly_n50(x))

    df = df.merge(df1.rename(columns={"assm": "source", "len": "source_len"}))
    df = df.merge(df1.rename(columns={"assm": "target", "len": "target_len"}))
    os.remove(rfile)
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

        log.debug(
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

            log.debug(
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

                log.debug(
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

            logging.debug(
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
                return x, weights[mid]
            else:

                log.debug(
                    "Couldn't find a filtering threshold, returning original graph"
                )
                return g, None

        if nx.is_connected(x):
            x, w = binary_search_filter(g, low=mid, high=high, weights=weights)
        else:
            x, w = binary_search_filter(g, low=low, high=mid, weights=weights)

        if nx.number_connected_components(x) == 1:
            return x, w
        else:
            log.debug("Couldn't find a filtering threshold, returning original graph")
            return g, None


def get_reps(graph, partition):
    log.debug("Finding the most central assemblies in the graph")
    """
    This function finds the representative genomes based on centrality measures
    """
    subgraphs = []
    reps = []
    for i in set([partition[k] for k in partition]):
        nodes = [k for k in partition if partition[k] == i]
        subgraph = graph.subgraph(nodes)
        subgraphs.append(subgraph)
        cents = nx.eigenvector_centrality(
            subgraph, weight="weight", max_iter=10000, tol=1e-4
        )
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


def create_slurm_commands(files, wdir, frag_len):
    file_list = list(product(files.keys(), files.keys()))
    odir = os.path.join(wdir, "fastANI_out")
    Path(odir).mkdir(parents=True, exist_ok=True)
    cmds = []
    ofiles = []
    for file in file_list:
        ofile = os.path.join(odir, "fastANI_out_" + str(file[0]) + "_" + str(file[1]))
        cmd = " ".join(
            [
                "fastANI",
                "--ql",
                files[file[0]],
                "--rl",
                files[file[1]],
                "--fragLen",
                str(int(frag_len)),
                "-o",
                ofile,
            ]
        )
        cmds.append(cmd)
        ofiles.append(ofile)
    return cmds, ofiles, odir


def check_slurm_output(files):
    pass


def split_fixed_size(lst, n):
    splits = [lst[i * n : (i + 1) * n] for i in range((len(lst) + n - 1) // n)]
    return splits


def map_slurm_jobs(cmds, slurm_config, odir, max_jobs_array, tmp_dir, threads):
    rfile = os.path.join(odir, "fastANI_out_jobs.txt")
    with open(rfile, "w") as outfile:
        for cmd in cmds:
            outfile.write("%s\n" % cmd)

    slurm = Slurm(
        **yaml.load(open(slurm_config), Loader=yaml.FullLoader),
        output=(
            # "{}/gderep-{}_{}.out".format(
            #     tmp_dir, Slurm.JOB_ARRAY_MASTER_ID, Slurm.JOB_ARRAY_ID
            "/dev/null"
            # )
        )
    )

    if len(cmds) > max_jobs_array:
        log.debug(
            "The number of jobs ({}) is larger than the allowed array size ({}). Splitting jobs..".format(
                len(cmds), max_jobs_array
            )
        )

    job_ranges = split_fixed_size(range(1, len(cmds) + 1), max_jobs_array)
    job_ids = []
    log.debug(
        "Mapping {} jobs to SLURM in {} batch(es)".format(len(cmds), len(job_ranges))
    )

    df = pd.read_csv(rfile, header=None)
    pw = []
    for r in job_ranges:
        rfile_tmp = os.path.join(odir, "fastANI_out_jobs-tmp.txt")
        skrows = min(r) - 1
        nrows = max(r)
        df1 = df.iloc[(min(r) - 1) : (max(r))]
        df1.columns = ["cmd"]

        ofiles = df1["cmd"].str.rsplit(" ", 1).str[-1].values

        df1.to_csv(rfile_tmp, index=False, sep="\t", header=False)
        new_r = range(1, df1.shape[0] + 1)

        slurm.set_array(new_r)
        slurm.add_arguments(wait="")
        bjob = "$(awk -vJ=" + "$SLURM_ARRAY_TASK_ID " + "'NR==J' " + rfile_tmp + " )"
        text_trap = io.StringIO()
        sys.stdout = text_trap
        job_id = slurm.sbatch(bjob)
        sys.stdout = sys.__stdout__
        job_ids.append(job_id)
        pairwise_distances = reduce_slurm_jobs(ofiles, threads)
        pw.append(pairwise_distances)
    return pw


def check_pw(pw, assms):
    res = {}

    if pw.empty:
        res = {"missing": [], "failed": True}
        return res

    df_source = pw.source.drop_duplicates().to_frame()
    df_source.columns = ["assm"]
    df_target = pw.target.drop_duplicates().to_frame()
    df_target.columns = ["assm"]
    df_assms = concat_df([df_source, df_target]).drop_duplicates()
    ids = df_assms["assm"].to_list()
    #    ids = list(set(pw["source"].tolist() + pw["target"].tolist()))
    # ids = list(set(pw.source).union(set(pw.target)))
    # Find if we are missing any assmebly (too divergent)
    diffs = list(set(assms) - set(ids))
    res = {"missing": diffs, "failed": False}
    return res


# From https://gist.github.com/TariqAHassan/fc77c00efef4897241f49e61ddbede9e
def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


def concat_df(frames):
    COLUMN_NAMES = frames[0].columns
    df_dict = dict.fromkeys(COLUMN_NAMES, [])
    for col in COLUMN_NAMES:
        extracted = (frame[col] for frame in frames)
        # Flatten and save to df_dict
        df_dict[col] = fast_flatten(extracted)
    df = pd.DataFrame.from_dict(df_dict)[COLUMN_NAMES]
    return df


def reduce_slurm_jobs(ofiles, threads):
    # if is_debug():
    #     dfs = list(
    #         map(process_fastANI_results, ofiles),
    #     )
    # else:
    p = Pool(threads)
    dfs = list(
        p.imap_unordered(process_fastANI_results, ofiles),
    )
    dfs = concat_df(dfs)
    return dfs


def is_unique(s):
    a = s.to_numpy()  # s.values (pandas<0.24)
    return (a[0] == a).all()


def estimate_frag_len(all_assemblies):
    assm_lens = (
        all_assemblies["assembly"].map(lambda x: get_assembly_length(x)).tolist()
    )
    x = min(assm_lens)
    if x < 3000:
        frag_len = math.floor(x / 500.0) * 500.0
    else:
        frag_len = 3000

    log.debug(
        "Minimum assembly length of {}. fastANI fragment length used: {}".format(
            x, frag_len
        )
    )

    if frag_len == 0:
        return None

    return frag_len


def create_graph(pairwise_distances):
    M = nx.from_pandas_edgelist(
        pairwise_distances, edge_attr=True, create_using=nx.MultiGraph()
    )
    G = nx.Graph()
    for u, v, data in M.edges(data=True):
        if not G.has_edge(u, v):
            # set weight to 1 if no weight is given for edge in M
            weight_raw = mean(
                d.get("weight_raw", 0) for d in M.get_edge_data(u, v).values()
            )
            weight = mean(d.get("weight", 0) for d in M.get_edge_data(u, v).values())
            G.add_edge(u, v, weight=weight, weight_raw=weight_raw)
    # Remove self loops
    G.remove_edges_from(list(nx.selfloop_edges(G)))
    # Identify isolated nodes
    isolated = list(nx.isolates(G))
    G.remove_nodes_from(isolated)
    return G.to_undirected(), isolated


def filter_graph(G):
    graphs_filt = []
    w_filts = []
    if G.number_of_edges() > 2:
        # TODO: Check how many components do we have
        n_comp = nx.number_connected_components(G)
        if n_comp > 1:
            log.debug("Graph with {} component(s".format(n_comp))
            for component in sorted(nx.connected_components(G), key=len, reverse=True):
                component = G.subgraph(component).copy()
                weights = []

                for u, v, d in component.edges(data=True):
                    weights.append(d["weight"])

                weights = sorted(set(weights))
                log.debug("Filtering graph")
                G_filt, w_filt = binary_search_filter(
                    g=component, low=0, high=len(weights) - 1, weights=weights
                )
                graphs_filt.append(G_filt)
                w_filts.append(w_filt)
        else:
            weights = []
            for u, v, d in G.edges(data=True):
                weights.append(d["weight"])
            weights = sorted(set(weights))
            log.debug("Filtering ANI graph")
            G_filt, w_filt = binary_search_filter(
                g=G, low=0, high=len(weights) - 1, weights=weights
            )
            graphs_filt.append(G_filt)
            w_filts.append(w_filt)
    else:
        log.debug("Skipping graph filtering, number of edges <=2")
        graphs_filt = [G]
        w_filts = [None]
    # partition = community_louvain.best_partition(G_filt, resolution=1.0)
    G_filt = nx.compose_all(graphs_filt)
    w_filt = ",".join(map(lambda x: str(x), w_filts))
    return G_filt, w_filt


def get_communities(G_filt):
    log.debug("Finding genome communities using Leiden community detection algorithm")
    g = __from_nx_to_igraph(G_filt, directed=False)
    part = leidenalg.find_partition(
        g,
        leidenalg.ModularityVertexPartition,
        weights="weight",
        seed=1234,
        n_iterations=-1,
    )
    leiden_coms = [g.vs[x]["name"] for x in part]
    partition = {k: i for i, sub_l in enumerate(leiden_coms) for k in sub_l}

    return partition


def create_mash_sketch(fname, temp_dir, threads):
    sketches = os.path.join(temp_dir, "mash.msh")

    mash_command = [
        "mash",
        "sketch",
        "-p",
        str(threads),
        "-o",
        temp_dir + "/mash",
        "-s",
        "10000",
        "-l",
        fname,
    ]
    log.debug("Generating MASH sketches")
    subprocess.run(mash_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return sketches


def run_mash(mash_sketch, threads, mash_threshold, temp_dir):

    mash_fout = os.path.join(temp_dir, "msh.out")
    mash_command = [
        "mash",
        "dist",
        "-p",
        str(threads),
        mash_sketch,
        mash_sketch,
    ]
    # mash_out = subprocess.run(mash_command, stdout=subprocess.PIPE).stdout.decode()
    log.debug("Running MASH")
    with open(mash_fout, "w") as fname:
        subprocess.run(mash_command, stdout=fname, stderr=subprocess.STDOUT)

    log.debug("Reading MASH results")
    mash_out = pd.read_csv(mash_fout, sep="\t", usecols=range(0, 3), header=None)
    mash_out.columns = ["source", "target", "weight"]
    mash_out.loc[:, "weight"] = 1 - mash_out["weight"]
    log.debug("Extracting comparisons with MASH distance >= {}".format(mash_threshold))
    mash_out = mash_out[mash_out["weight"] >= mash_threshold]
    log.debug("Calculating genome lengths")

    # df = pd.DataFrame(mash_out, columns=["source", "target", "weight"])
    # df["weight"] = df["weight"].astype(float)
    # df["weight"] = df["weight"].apply(lambda x: 1.0 - x)

    df_source = mash_out.source.drop_duplicates().reset_index(drop=True)
    df_target = mash_out.target.drop_duplicates().reset_index(drop=True)
    df_assms = pd.concat([df_source, df_target]).to_frame()
    df_assms.columns = ["assm"]
    p = ThreadPool(processes=threads)
    df_assms.loc[:, "len"] = list(
        p.imap(get_assembly_length, df_assms["assm"].to_list())
    )

    df_assms.columns = ["source", "source_len"]
    mash_out = mash_out.merge(df_assms)

    df_assms.columns = ["target", "target_len"]
    mash_out = mash_out.merge(df_assms)

    return mash_out


def dereplicate_mash(all_assemblies, threads, tmp_dir, mash_threshold, threshold):
    partitions = []
    all_assemblies_tmp = all_assemblies.copy()
    w_filt = 1 - mash_threshold

    with tempfile.TemporaryDirectory(dir=tmp_dir, prefix="gderep-mash") as temp_dir:
        fname = os.path.join(temp_dir, "assemblies.txt")
        all_assemblies_tmp["assembly"].to_csv(fname, index=False, header=False)
        mash_sketch = create_mash_sketch(fname, temp_dir, threads)
        mash_dist = run_mash(
            mash_sketch, threads=threads, mash_threshold=w_filt, temp_dir=temp_dir
        )

    log.debug("Checking we have all genomes")
    pw = check_pw(pw=mash_dist, assms=all_assemblies_tmp["assembly"].tolist())

    if pw["failed"]:
        log.debug("Failed. Reason: Assemblies too divergent")
        failed = "Assemblies too divergent"
        return None, None, failed

    missing = pw["missing"]

    if len(missing) > 0:
        log.debug(
            "Missing {:,} assemblies in the MASH comparison, Assemblies too divergent".format(
                len(missing)
            )
        )

    # Convert pw dist to graph
    log.debug("Generating MASH graph")
    # mash_dist = mash_dist[(mash_dist.weight >= float(w_filt))]

    G, isolated = create_graph(mash_dist)
    if G.number_of_edges() == 0:
        log.debug("Only self loops found")
        all_assemblies_tmp.loc[:, "representative"] = 1
        all_assemblies_tmp.loc[:, "derep"] = 1

        return all_assemblies_tmp

    # TODO: streamline this
    # G_filt, w_filt = filter_graph(G)
    # edges2rm = [
    #     (u, v) for u, v, d in G.edges(data=True) if float(d["weight"]) <= float(w_filt)
    # ]
    # G.remove_edges_from(edges2rm)
    G_filt = G
    partition = get_communities(G_filt)

    log.debug(
        "Graph properties: w_filt={} edges={} density={} components={} communities={}".format(
            w_filt,
            G_filt.number_of_edges(),
            nx.density(G_filt),
            nx.number_connected_components(G_filt),
            len(set([partition[k] for k in partition])),
        )
    )
    reps, subgraphs = get_reps(graph=G_filt, partition=partition)
    derep_assemblies = []
    for i in range(len(reps)):
        # candidates = refine_candidates(
        #     rep=reps[i],
        #     subgraph=subgraphs[i],
        #     threshold=threshold,
        #     pw=mash_dist,
        # )
        can = representative(
            name=reps[i], graph=subgraphs[i], threshold=threshold, pw=mash_dist
        )

        candidates = [can.name] + can.selected

        if len(candidates) > 1:
            for candidate in candidates:
                derep_assemblies.append(candidate)
        else:
            derep_assemblies = candidates

    derep_assemblies = derep_assemblies + reps + missing + isolated
    derep_assemblies = list(set(derep_assemblies))
    reps = reps + missing + isolated

    if len(reps) > len(derep_assemblies):
        log.debug(
            "Less representatives {} than dereplicated assemblies {} !!!".format(
                len(reps), len(derep_assemblies)
            )
        )
        exit(0)

    # all_assemblies = all_assemblies[
    #     all_assemblies["assembly"].isin(derep_assemblies)
    # ]
    all_assemblies_tmp.loc[:, "representative"] = 0
    all_assemblies_tmp.loc[
        all_assemblies_tmp["assembly"].isin(reps), "representative"
    ] = 1
    all_assemblies_tmp.loc[:, "derep"] = 0
    all_assemblies_tmp.loc[
        all_assemblies_tmp["assembly"].isin(derep_assemblies), "derep"
    ] = 1
    return all_assemblies_tmp


def dereplicate_ANI(
    all_assemblies, threads, threshold, chunks, slurm_config, tmp_dir, max_jobs_array
):
    """
    This function dereplicates genomes by Taxon by:
    1. Calculating ANI and fraction aligned
    2. Dynamically filters the ANI graph
    3.
    """

    derep_assemblies = []
    n_assemblies = all_assemblies.shape[0]
    all_assemblies_tmp = all_assemblies.copy()
    with tempfile.TemporaryDirectory(dir=tmp_dir, prefix="gderep-") as temp_dir:
        failed = None
        if n_assemblies < 10 or slurm_config is None:
            if (n_assemblies * n_assemblies) < threads:
                threads = n_assemblies * n_assemblies
            log.debug(
                "Found {} assemblies, using default fastANI with {} threads".format(
                    n_assemblies, threads
                )
            )

            frag_len = estimate_frag_len(all_assemblies_tmp)

            if frag_len is None:
                log.debug("Failed. Reason: Assemblies too short")
                failed = "Assemblies too short"
                return None, None, failed, None

            ani_results = pairwise_fastANI(
                assemblies=all_assemblies_tmp["assembly"].tolist(),
                threads=threads,
                temp_dir=temp_dir,
                frag_len=frag_len,
            )

            log.debug("Processing ANI results")
            pairwise_distances = process_fastANI_results(ani_results)
            log.debug("Obtained {:,} comparisons".format(pairwise_distances.shape[0]))
        else:
            n = chunks
            chunks = split_fixed_size(all_assemblies_tmp["assembly"].tolist(), n)

            log.debug(
                "Found {} assemblies, creating {} chunks with {} genomes".format(
                    len(all_assemblies_tmp), len(chunks), n
                )
            )

            frag_len = estimate_frag_len(all_assemblies_tmp)

            if frag_len is None:
                log.debug("Failed. Reason: Assemblies too short")
                failed = "Assemblies too short"
                return None, None, failed, None

            # save chunks to disk
            files, wdir = save_chunks_to_disk(chunks, temp_dir)
            # Create fastANI commands
            cmds, ofiles, odir = create_slurm_commands(files, wdir, frag_len)
            # Run slurm array job
            pairwise_distances = map_slurm_jobs(
                cmds, slurm_config, odir, max_jobs_array, temp_dir, threads
            )
            pairwise_distances = concat_df(pairwise_distances)
            log.debug("Reducing SLURM jobs and processing ANI results")
            # pairwise_distances = reduce_slurm_jobs(ofiles, threads)

        pw = check_pw(
            pw=pairwise_distances, assms=all_assemblies_tmp["assembly"].tolist()
        )

        if pw["failed"]:
            log.debug("Failed. Reason: Assemblies too divergent")
            failed = "Assemblies too divergent"
            return None, None, failed, None

        missing = pw["missing"]

        if len(missing) > 0:
            log.debug(
                "Missing {:,} assemblies in the fastANI comparison, Assemblies too divergent".format(
                    len(missing)
                )
            )

        # pairwise_distances[
        #     ["source", "target", "ANI", "source_len", "target_len"]
        # ].to_csv("test.tsv", index=False, sep="\t")
        # Convert pw dist to graph
        log.debug("Generating ANI graph")

        G, isolated = create_graph(pairwise_distances)

        if G.number_of_edges() == 0:
            log.debug("Only self loops found")
            all_assemblies_tmp.loc[:, "representative"] = 1
            all_assemblies_tmp.loc[:, "derep"] = 1

            results = pd.DataFrame(
                [
                    (
                        0,
                        G.number_of_nodes(),
                        G.number_of_nodes(),
                        G.number_of_edges(),
                    )
                ],
                columns=["weight", "communities", "n_genomes", "n_genomes_derep"],
            )

            return all_assemblies_tmp, results, failed, None

        # TODO: streamline this
        G_filt, w_filt = filter_graph(G)

        partition = get_communities(G_filt)

        log.debug(
            "Graph properties: w_filt={} edges={} density={} components={} communities={}".format(
                w_filt,
                G_filt.number_of_edges(),
                nx.density(G_filt),
                nx.number_connected_components(G_filt),
                len(set([partition[k] for k in partition])),
            )
        )
        log.debug(
            "Refining genome selection (length difference and fraction aligned, z-score={})".format(
                threshold
            )
        )

        reps, subgraphs = get_reps(graph=G_filt, partition=partition)
        derep_assemblies = []
        stats = []
        for i in range(len(reps)):
            # candidates = refine_candidates(
            #     rep=reps[i],
            #     subgraph=subgraphs[i],
            #     threshold=threshold,
            #     pw=pairwise_distances,
            # )
            can = representative(
                name=reps[i],
                graph=subgraphs[i],
                threshold=threshold,
                pw=pairwise_distances,
            )
            stats_df = pd.DataFrame(
                [
                    (
                        reps[i],
                        can.n_nodes,
                        can.n_nodes_selected,
                        can.n_nodes_discarded,
                        can.graph_avg_weight,
                        can.graph_sd_weight,
                        can.graph_avg_weight_raw,
                        can.graph_sd_weight_raw,
                        can.subgraph_selected_avg_weight,
                        can.subgraph_selected_sd_weight,
                        can.subgraph_selected_avg_weight_raw,
                        can.subgraph_selected_sd_weight_raw,
                        can.subgraph_discarded_avg_weight,
                        can.subgraph_discarded_sd_weight,
                        can.subgraph_discarded_avg_weight_raw,
                        can.subgraph_discarded_sd_weight_raw,
                    )
                ],
                columns=[
                    "representative",
                    "n_nodes",
                    "n_nodes_selected",
                    "n_nodes_discarded",
                    "graph_avg_weight",
                    "graph_sd_weight",
                    "graph_avg_weight_raw",
                    "graph_sd_weight_raw",
                    "subgraph_selected_avg_weight",
                    "subgraph_selected_sd_weight",
                    "subgraph_selected_avg_weight_raw",
                    "subgraph_selected_sd_weight_raw",
                    "subgraph_discarded_avg_weight",
                    "subgraph_discarded_sd_weight",
                    "subgraph_discarded_avg_weight_raw",
                    "subgraph_discarded_sd_weight_raw",
                ],
            )

            stats.append(stats_df)
            candidates = [can.name] + can.selected

            if len(candidates) > 1:
                for candidate in candidates:
                    derep_assemblies.append(candidate)
            else:
                derep_assemblies = candidates
    stats = pd.concat(stats)
    # derep_assemblies = [item for derep_assemblies in t for item in derep_assemblies]
    derep_assemblies = derep_assemblies + reps + missing + isolated

    derep_assemblies = list(set(derep_assemblies))
    reps = reps + missing + isolated

    log.debug(
        "Keeping {}/{} genomes".format(len(derep_assemblies), len(all_assemblies_tmp))
    )
    results = pd.DataFrame(
        [
            (
                w_filt,
                len(set([partition[k] for k in partition])),
                len(all_assemblies_tmp),
                len(derep_assemblies),
            )
        ],
        columns=["weight", "communities", "n_genomes", "n_genomes_derep"],
    )

    if len(reps) > len(derep_assemblies):
        log.debug(
            "Less representatives {} than dereplicated assemblies {} !!!".format(
                len(reps), len(derep_assemblies)
            )
        )
        exit(0)

    rep_keys = all_assemblies_tmp[all_assemblies_tmp["assembly"].isin(reps)][
        "accession"
    ].tolist()

    # all_assemblies = all_assemblies[
    #     all_assemblies["assembly"].isin(derep_assemblies)
    # ]
    all_assemblies_tmp.loc[:, "representative"] = 0
    all_assemblies_tmp.loc[
        all_assemblies_tmp["assembly"].isin(reps), "representative"
    ] = 1
    all_assemblies_tmp.loc[:, "derep"] = 0
    all_assemblies_tmp.loc[
        all_assemblies_tmp["assembly"].isin(derep_assemblies), "derep"
    ] = 1

    return all_assemblies_tmp, results, failed, stats


# def refine_candidates(rep, subgraph, pw, threshold=2.0):
#     results = []
#     len_diff = 0
#     source_len = 0
#     if subgraph.number_of_edges() == 0:
#         return [rep]
#     for reachable_node in nx.all_neighbors(subgraph, rep):
#         if reachable_node != rep:
#             idx = pw["target"] == rep
#             pw.loc[idx, ["source", "target", "source_len", "target_len"]] = pw.loc[
#                 idx, ["target", "source", "target_len", "source_len"]
#             ].values

#             df = pw[(pw["source"] == rep) & (pw["target"] == reachable_node)]

#             # Switch column order
#             # idx = df["target"] == rep
#             # df.loc[idx, ["source", "target", "source_len", "target_len"]] = df.loc[
#             #     idx, ["target", "source", "target_len", "source_len"]
#             # ].values

#             source_len = mean(df.source_len.values)
#             target_len = mean(df.target_len.values)
#             len_diff = abs(source_len - target_len)
#             weight = mean(df.weight.values)

#             if "aln_frac" in df.columns:
#                 aln_frac = mean(df.aln_frac.values)
#                 results.append(
#                     {
#                         "source": rep,
#                         "target": reachable_node,
#                         "weight": float(weight),
#                         "aln_frac": float(aln_frac),
#                         "len_diff": int(len_diff),
#                     }
#                 )
#             else:
#                 results.append(
#                     {
#                         "source": rep,
#                         "target": reachable_node,
#                         "weight": float(weight),
#                         "len_diff": int(len_diff),
#                     }
#                 )
#     df = pd.DataFrame(results)

#     if len(df.index) < 2:
#         if len_diff > 0:
#             diff_ratio = source_len / len_diff
#         else:
#             diff_ratio = 0

#         if diff_ratio >= 0.1:
#             assms = df["target"].tolist()
#             assms.append(rep)
#         else:
#             assms = [rep]
#     else:
#         if "aln_frac" in df.columns:
#             if is_unique(df["aln_frac"]):
#                 df["z_score_aln"] = 0
#             else:
#                 df["z_score_aln"] = stats.zscore(df["aln_frac"])

#         if is_unique(df["len_diff"]):
#             df["z_score_diff"] = 0
#         else:
#             df["z_score_diff"] = stats.zscore(df["len_diff"])

#         if "aln_frac" in df.columns:
#             df = df[
#                 (df["z_score_aln"].abs() >= threshold)
#                 | (df["z_score_diff"].abs() >= threshold)
#             ]
#         else:
#             df = df[(df["z_score_diff"].abs() >= threshold)]

#         assms = df["target"].tolist()
#         assms.append(rep)

#     if df.empty:
#         return [rep]
#     else:
#         return assms
