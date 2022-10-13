"""
This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

import subprocess
import sys
import os
import shutil
from derep_genomes.general import (
    get_arguments,
    find_assemblies_for_accessions,
    is_debug,
)
from derep_genomes.graph import dereplicate_ANI, dereplicate_xash, concat_df
import logging
import pathlib
import sqlite3
from derep_genomes.dbops import (
    check_if_db_exists,
    check_if_db_empty,
    create_db_tables,
    check_db_tables,
    check_if_done,
    retrieve_all_taxa_analyzed,
    retrieve_all_jobs_done,
    retrieve_all_jobs_failed,
    delete_from_db,
    retrieve_all_genomes_derep,
)
from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor, as_completed

from tqdm import tqdm
from functools import partial
import pandas as pd

log = logging.getLogger("my_logger")


def process_one_taxon(taxon, parms):
    classification = parms["classification"]
    threads = parms["threads"]
    threshold = parms["threshold"]
    chunks = parms["chunks"]
    slurm_config = parms["slurm_config"]
    slurm_threads = parms["slurm_threads"]
    tmp_dir = parms["tmp_dir"]
    max_jobs_array = parms["max_jobs_array"]
    xash_threshold = parms["xash_threshold"]
    ani_fraglen_fraction = parms["ani_fraglen_fraction"]
    min_genome_size = parms["min_genome_size"]
    failed_df = None
    dashing = parms["dashing"]
    assm_max = parms["assm_max"]
    threads = parms["threads"]

    log.debug("Dereplicating {}".format(taxon))
    classification = classification[classification["taxonomy"] == taxon]
    acc_to_assemblies = classification.loc[:, ["accession", "file"]]
    n_assemblies = classification.shape[0]
    # check if already processed

    if n_assemblies == 0:
        return
    else:
        log.debug(
            "{:,} assemblies for this species, clustering to dereplicate.".format(
                n_assemblies
            )
        )
        if dashing:
            log.debug("Preclustering with Dashing...")
        else:
            log.debug("Preclustering with Mash...")

        acc_to_assemblies_xash = dereplicate_xash(
            all_assemblies=acc_to_assemblies,
            threads=threads,
            tmp_dir=tmp_dir,
            xash_threshold=xash_threshold,
            threshold=threshold,
            dashing=dashing,
        )
        acc_to_assemblies_xash_r = acc_to_assemblies_xash[
            acc_to_assemblies_xash["representative"] == 1
        ].copy()
        acc_to_assemblies_xash_d = acc_to_assemblies_xash[
            acc_to_assemblies_xash["derep"] == 1
        ].copy()
        acc_to_assemblies_xash_o = acc_to_assemblies_xash[
            acc_to_assemblies_xash["derep"] == 0
        ].copy()

        if dashing:
            log.debug(
                "Identified {:,} linked assemblies in {:,} Dashing communities".format(
                    acc_to_assemblies_xash_d.shape[0], acc_to_assemblies_xash_r.shape[0]
                )
            )
        else:
            log.debug(
                "Identified {:,} linked assemblies in {:,} MASH communities".format(
                    acc_to_assemblies_xash_d.shape[0], acc_to_assemblies_xash_r.shape[0]
                )
            )

        # We refine each MASH community
        if dashing:
            log.debug(
                "Refining {:,} Dashing selected assemblies".format(
                    acc_to_assemblies_xash_d.shape[0]
                )
            )
        else:
            log.debug(
                "Refining {:,} MASH communities".format(
                    acc_to_assemblies_xash_r.shape[0]
                )
            )
        derep_assemblies, results, failed, stats_df = dereplicate_ANI(
            all_assemblies=acc_to_assemblies_xash_d,
            threads=threads,
            threshold=threshold,
            chunks=chunks,
            slurm_config=slurm_config,
            slurm_threads=slurm_threads,
            tmp_dir=tmp_dir,
            max_jobs_array=max_jobs_array,
            min_genome_size=min_genome_size,
            ani_fraglen_fraction=ani_fraglen_fraction,
            assm_max=assm_max,
        )

        if failed is not None:
            failed_df = acc_to_assemblies.copy()
            failed_df.loc[:, "file"] = failed_df["file"]
            failed_df.loc[:, "taxonomy"] = taxon
            failed_df.loc[:, "reason"] = failed
            failed_df = failed_df[["taxonomy", "accession", "file", "reason"]]
        else:
            failed_df = pd.DataFrame()
            # Add taxa to taxa table
            derep_assemblies = concat_df([derep_assemblies, acc_to_assemblies_xash_o])
            derep_assemblies.loc[:, "taxonomy"] = taxon
            results.loc[:, "taxonomy"] = taxon
            if stats_df is not None:
                stats_df.loc[:, "taxonomy"] = taxon
            else:
                stats_df = pd.DataFrame()

        return derep_assemblies, results, failed_df, stats_df


def insert_to_db(
    derep_assemblies, results, failed, stats_df, con, threads, out_dir, copy
):
    if not derep_assemblies.empty:
        query = check_existing(
            df=derep_assemblies.loc[:, ["taxonomy"]].drop_duplicates(),
            columns=["taxonomy"],
            table="taxa",
            con=con,
        )
        insert_pd_sql(query=query, table="taxa", con=con)

        # Add genomes to genome table
        query = check_existing(
            df=derep_assemblies.loc[:, ["taxonomy", "accession"]],
            columns=["taxonomy", "accession"],
            table="genomes",
            con=con,
        )
        genomes = insert_pd_sql(query=query, table="genomes", con=con)

        # Add derep genomes to derep genome table
        query = check_existing(
            df=derep_assemblies[derep_assemblies["derep"] == 1]
            .loc[:, ["accession", "representative"]]
            .copy(),
            columns=["accession", "representative"],
            table="genomes_derep",
            con=con,
        )
        genomes_derep = insert_pd_sql(query=query, table="genomes_derep", con=con)

    if not results.empty:
        # Add results to results table
        query = check_existing(
            results.loc[
                :, ["taxonomy", "weight", "communities", "n_genomes", "n_genomes_derep"]
            ],
            columns=[
                "taxonomy",
                "weight",
                "communities",
                "n_genomes",
                "n_genomes_derep",
            ],
            table="results",
            con=con,
        )
        results_db = insert_pd_sql(query=query, table="results", con=con)

    if not stats_df.empty:
        # Add results to results table
        query = check_existing(
            stats_df.loc[
                :,
                [
                    "taxonomy",
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
            ],
            columns=[
                "taxonomy",
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
            table="stats",
            con=con,
        )

        stats_db = insert_pd_sql(query=query, table="stats", con=con)

    if not derep_assemblies.empty and not results.empty:
        # Add jobs done to table
        query = check_existing(
            derep_assemblies.loc[:, ["taxonomy", "accession", "file"]],
            columns=["taxonomy", "accession", "file"],
            table="jobs_done",
            con=con,
        )
        jobs_done_db = insert_pd_sql(query=query, table="jobs_done", con=con)

    # Add failed jobs to table
    if not failed.empty:
        query = check_existing(
            failed.loc[:, ["taxonomy", "accession", "file", "reason"]],
            columns=["taxonomy", "accession", "file", "reason"],
            table="jobs_failed",
            con=con,
        )
        insert_pd_sql(query=query, table="jobs_failed", con=con)

        query = check_existing(
            df=failed.loc[:, ["taxonomy"]].drop_duplicates(),
            columns=["taxonomy"],
            table="taxa",
            con=con,
        )
        insert_pd_sql(query=query, table="taxa", con=con)

        query = check_existing(
            failed.loc[:, ["taxonomy", "accession", "file"]],
            columns=["taxonomy", "accession", "file"],
            table="jobs_done",
            con=con,
        )
        insert_pd_sql(query=query, table="jobs_done", con=con)

    try:
        con.commit()
    except:
        con.rollback()


def copy_files(src_dst):
    log.debug("{} -> {}".format(src_dst[0], src_dst[1]))
    shutil.copy(src_dst[0], src_dst[1])


def symlink_files(src, dst):
    os.symlink(src, os.path.join(dst, os.path.basename(src)))


def check_existing(df, columns, table, con):
    query = "SELECT " + ",".join(columns) + " from " + table
    old = pd.read_sql_query(query, con)
    new = df.loc[:, columns]
    old = old.loc[:, columns]
    if old.empty:
        return new
    else:
        return new.loc[~new.isin(old).all(1)]


def insert_pd_sql(query, table, con):
    if query.empty:
        log.debug("All entries in table {} already present".format(table))
    else:
        log.debug("Adding {} entries in table {}".format(query.shape[0], table))
        query.to_sql(name=table, con=con, if_exists="append", index=False)


def process_sigletons(singletons, out_dir, threads, con, copy):

    # Insert data to table taxa
    query = check_existing(singletons, columns=["taxonomy"], table="taxa", con=con)
    insert_pd_sql(query=query, table="taxa", con=con)

    # Insert data to table genomes
    query = check_existing(
        singletons, columns=["taxonomy", "accession"], table="genomes", con=con
    )
    genomes = insert_pd_sql(query=query, table="genomes", con=con)

    # Insert data to table genomes_derep
    singletons.loc[:, "representative"] = int(1)

    query = check_existing(
        singletons,
        columns=["accession", "representative"],
        table="genomes_derep",
        con=con,
    )
    genomes_derep = insert_pd_sql(query=query, table="genomes_derep", con=con)

    # Insert data to table results
    singletons.loc[:, "weight"] = None
    singletons.loc[:, "communities"] = int(1)
    singletons.loc[:, "n_genomes"] = 1
    singletons.loc[:, "n_genomes_derep"] = 1

    query = check_existing(
        singletons,
        columns=["taxonomy", "weight", "communities", "n_genomes", "n_genomes_derep"],
        table="results",
        con=con,
    )
    results = insert_pd_sql(query=query, table="results", con=con)
    # Insert data to table jobs_done
    # singletons.loc[:, "file"] = singletons["assembly"].apply(
    #     lambda x: os.path.basename(x)
    # )
    query = check_existing(
        singletons,
        columns=["taxonomy", "accession", "file"],
        table="jobs_done",
        con=con,
    )
    jobs_done = insert_pd_sql(query=query, table="jobs_done", con=con)

    try:
        con.commit()
    except:
        pass


def find_assemblies(x, classifications, all_assm):
    accessions = sorted(classifications[x])
    log.debug(x)
    log.debug(accessions[0])
    res = find_assemblies_for_accessions(accessions=accessions, all_assemblies=all_assm)
    return res


# TODO: Improve the matching between filenames and accessions
def shorten_accession(accession):
    if accession.startswith("GCF_") or accession.startswith("GCA_"):
        accession = accession.split(".")[0]
        assert len(accession) == 13
    else:
        accession = accession
    return accession


def get_accession(assembly):
    res = {}
    accession = os.path.basename(assembly)
    accession = shorten_accession(accession)
    res = {"accession_nover": accession, "assembly": assembly}
    return res


def check_done(data, con):
    df = check_existing(
        data, columns=["taxonomy", "accession", "file"], table="jobs_done", con=con
    )
    return df


def check_done_apply(df, parms):
    db = parms["db"]
    con = sqlite3.connect(db, isolation_level="EXCLUSIVE")
    taxon = list(set(df["taxonomy"]))[0]
    acc_to_assemblies = df.loc[:, ["accession", "file"]]
    is_done = check_if_done(con=con, taxon=taxon, acc2assm=acc_to_assemblies)

    con.close()
    return is_done


def command_exists(command):
    try:
        fnull = open(os.devnull, "w")
        subprocess.call([command], stdout=fnull, stderr=subprocess.STDOUT)
        return True
    except OSError:
        return False


def mute():
    sys.stdout = open(os.devnull, "w")


def file_exists(file):
    return os.path.isfile(file)


def main():

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    args = get_arguments()
    logging.getLogger("my_logger").setLevel(
        logging.DEBUG if args.debug else logging.INFO
    )
    prefix = args.prefix

    # test if binaries for fastani and mash exists
    # check if command exists
    if not command_exists("fastANI"):
        log.error("fastANI not found in PATH")
        sys.exit(1)
    if not command_exists("mash"):
        log.error("mash not found in PATH")
        sys.exit(1)
    if not command_exists("dashing"):
        log.error("dashing not found in PATH")
        sys.exit(1)
    if args.dashing:
        log.debug("Using dashing")
    else:
        log.debug("Using mash")

    if args.copy:
        out_dir = pathlib.Path(args.out_dir).absolute()
        os.makedirs(out_dir, exist_ok=True)

    tmp_dir = pathlib.Path(args.tmp_dir).absolute()
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir, exist_ok=True)

    con = check_if_db_exists(args.db)
    db_empty = check_if_db_empty(con)

    if db_empty:
        log.info("Creating db tables")
        create_db_tables(con)
    else:
        log.info("Checking correct tables exist")
        check_db_tables(con)

    # Get all assemblies
    log.info("Reading genome data")
    assm_data = pd.read_csv(args.tax_file, sep="\t")
    log.info("Found {:,} assemblies".format(assm_data.shape[0]))
    # Check if all files are present
    log.info("Checking if all files are present")

    p = Pool(processes=args.threads)
    file_list = assm_data["file"].tolist()

    file_list = list(
        tqdm(
            p.imap(os.path.isfile, file_list),
            total=len(file_list),
            leave=False,
            ncols=80,
            desc="Files processed",
        )
    )
    p.close()
    p.join()

    assm_data["file_exists"] = file_list
    p.close()
    p.join()

    # remove files that are not present
    pre_filt_assm = assm_data.shape[0]
    assm_data = assm_data[assm_data["file_exists"] == True]
    if pre_filt_assm != assm_data.shape[0]:
        log.info(
            "Removed {:,} assemblies that are not present".format(
                pre_filt_assm - assm_data.shape[0]
            )
        )
    else:
        log.info("All files are present")

    if args.selected_taxa:
        with args.selected_taxa as f:
            taxas = [line.rstrip() for line in f]
        assm_data = assm_data[assm_data["taxonomy"].isin(taxas)]

    taxa_in_db = retrieve_all_taxa_analyzed(con)
    jobs_done = retrieve_all_jobs_done(con)
    jobs_done_old = jobs_done.copy()

    log.info("Found {:,} taxa in the DB".format(taxa_in_db.shape[0]))

    log.info("Checking assemblies already processed")

    if not jobs_done.empty and not taxa_in_db.empty:
        # check that there are no updates
        jobs_done = jobs_done.merge(taxa_in_db).loc[:, ["accession", "file"]]
        to_do = pd.concat(
            [assm_data.loc[:, ["accession", "file"]], jobs_done, jobs_done]
        ).drop_duplicates(keep=False)

        to_do = assm_data[assm_data["accession"].isin(to_do["accession"])]

        log.info("{:,} assemblies already processed\n".format(jobs_done.shape[0]))
        to_rm = pd.concat([to_do.loc[:, ["taxonomy"]].drop_duplicates(), taxa_in_db])
        to_rm = to_rm[to_rm["taxonomy"].duplicated()].reset_index()

        if not to_rm.empty:
            log.info("Found {:,} taxa that need to be updated\n".format(to_rm.shape[0]))
            removed = delete_from_db(taxons=to_rm, con=con)
            log.warning(
                "Removed {:,} taxa and {:,} assemblies from DB".format(
                    to_rm.shape[0], removed.shape[0]
                )
            )
            if args.copy:
                files_to_rm = (
                    removed["file"].apply(lambda x: os.path.join(out_dir, x)).tolist()
                )
                if files_to_rm:
                    files_removed = []
                    for file in files_to_rm:
                        if os.path.isfile(file):
                            os.remove(file)
                            files_removed.append(file)
                    log.warning(
                        "Removed {} assemblies from {}\n".format(
                            len(files_removed), out_dir
                        )
                    )

            rfname = prefix + "-derep-genomes_removed-db-entries.tsv"
            removed.to_csv(rfname, index=False, sep="\t")

    else:
        log.info("No jobs found in the database\n")
        to_do = assm_data

    if not to_do.empty:
        assm_min = 2
        assm_max = args.assm_max

        log.info("Splitting taxa in between singleton and non-singleton taxa")

        taxon_counts = to_do.groupby("taxonomy", as_index=False)["taxonomy"].agg(
            {"count": "count"}
        )

        classification_singletons = taxon_counts[taxon_counts["count"] < 2][
            "taxonomy"
        ].tolist()
        classification_singletons = to_do[
            to_do["taxonomy"].isin(classification_singletons)
        ].reset_index(drop=True)

        if args.slurm_config:
            classification_small = taxon_counts[
                taxon_counts["count"].isin(range(assm_min, assm_max + 1))
            ]["taxonomy"].tolist()
            classification_small = to_do[
                to_do["taxonomy"].isin(classification_small)
            ].reset_index(drop=True)

            classification_large = taxon_counts[taxon_counts["count"] > assm_max][
                "taxonomy"
            ].tolist()

            classification_large = to_do[
                to_do["taxonomy"].isin(classification_large)
            ].reset_index(drop=True)
        else:
            classification_small = taxon_counts[taxon_counts["count"] > 1][
                "taxonomy"
            ].tolist()
            classification_small = to_do[
                to_do["taxonomy"].isin(classification_small)
            ].reset_index(drop=True)
            classification_large = pd.DataFrame()

        if not classification_singletons.empty:
            log.info(
                "Processing {:,} singletons".format(classification_singletons.shape[0])
            )
            process_sigletons(
                singletons=classification_singletons,
                out_dir=args.out_dir,
                threads=args.threads,
                con=con,
                copy=args.copy,
            )
        else:
            log.info("No singletons found\n")

        if not classification_small.empty:
            if args.workers > 2:
                workers = args.workers
                processes = args.threads // workers
            else:
                workers = 1
                processes = args.threads
            taxons = list(set(classification_small["taxonomy"]))
            log.info(
                "Dereplicating {:,} taxa with {} to {} assemblies using {} workers with {} threads".format(
                    len(taxons), assm_min, assm_max, workers, processes
                )
            )
            parms_small = {
                "classification": classification_small,
                "threads": processes,
                "threshold": args.threshold,
                "chunks": args.chunks,
                "slurm_config": None,
                "slurm_threads": None,
                "tmp_dir": tmp_dir,
                "max_jobs_array": args.max_jobs_array,
                "xash_threshold": args.xash_threshold,
                "min_genome_size": args.min_genome_size,
                "ani_fraglen_fraction": args.ani_fraglen_fraction,
                "dashing": args.dashing,
                "assm_max": assm_max,
                "out_dir": args.out_dir,
                # "con": con,
                "copy": args.copy,
            }

            if is_debug():
                dfs = list(
                    map(partial(process_one_taxon, parms=parms_small), taxons),
                )
            else:

                pbar = tqdm(
                    total=len(taxons),
                    desc="Taxa processed",
                    leave=False,
                    ncols=80,
                )  # Init pbar
                with ProcessPoolExecutor(workers) as executor:
                    futures = [
                        executor.submit(
                            partial(process_one_taxon, taxon=tax, parms=parms_small)
                        )
                        for tax in taxons
                    ]
                    for _ in as_completed(futures):
                        pbar.update(n=1)  # Increments counter
                dfs = [f.result() for f in futures]
                pbar.close()
            # Get all failed jobs
            failed_list = [x[2] for x in dfs if None not in dfs]
            if all(v is None for v in failed_list):
                failed = pd.DataFrame()
            else:
                failed = pd.concat([x[2] for x in dfs if None not in dfs])

            derep_assemblies_list = [x[0] for x in dfs if None not in dfs]
            if all(v is None for v in derep_assemblies_list):
                derep_assemblies = pd.DataFrame()
            else:
                derep_assemblies = pd.concat([x[0] for x in dfs if None not in dfs])

            results_list = [x[1] for x in dfs if None not in dfs]
            if all(v is None for v in results_list):
                results = pd.DataFrame()
            else:
                results = pd.concat([x[1] for x in dfs if None not in dfs])

            stats_df_list = [x[3] for x in dfs if None not in dfs]
            if all(v is None for v in stats_df_list):
                stats_df = pd.DataFrame()
            else:
                stats_df = pd.concat([x[3] for x in dfs if None not in dfs])

            # Insert results into database
            insert_to_db(
                derep_assemblies=derep_assemblies,
                results=results,
                failed=failed,
                stats_df=stats_df,
                con=con,
                threads=args.threads,
                out_dir=args.out_dir,
                copy=args.copy,
            )
        else:
            log.info("No taxa with 2 to {} assemblies found\n".format(assm_max))

        if not classification_large.empty:
            taxons = list(set(classification_large["taxonomy"]))
            log.info(
                "Dereplicating {:,} taxa with more than {} assemblies using SLURM with {} threads".format(
                    len(taxons), assm_max, args.slurm_threads
                )
            )
            parms_large = {
                "classification": classification_large,
                "threads": args.threads,
                "threshold": args.threshold,
                "chunks": args.chunks,
                "slurm_config": args.slurm_config.name,
                "slurm_threads": args.slurm_threads,
                "tmp_dir": tmp_dir,
                "max_jobs_array": args.max_jobs_array,
                "xash_threshold": args.xash_threshold,
                "min_genome_size": args.min_genome_size,
                "ani_fraglen_fraction": args.ani_fraglen_fraction,
                "dashing": args.dashing,
                "assm_max": assm_max,
                "out_dir": args.out_dir,
                "copy": args.copy,
            }
            if is_debug():
                dfs = list(
                    map(partial(process_one_taxon, parms=parms_large), taxons),
                )
            else:
                dfs = list(
                    tqdm(
                        map(partial(process_one_taxon, parms=parms_large), taxons),
                        total=len(taxons),
                        leave=False,
                        ncols=80,
                    )
                )
        else:
            if args.slurm_config:
                log.info(
                    "No taxa with more than {} assemblies found\n".format(assm_max)
                )

            # Get all failed jobs
            failed_list = [x[2] for x in dfs if None not in dfs]
            if all(v is None for v in failed_list):
                failed = pd.DataFrame()
            else:
                failed = pd.concat([x[2] for x in dfs if None not in dfs])

            derep_assemblies_list = [x[0] for x in dfs if None not in dfs]
            if all(v is None for v in derep_assemblies_list):
                derep_assemblies = pd.DataFrame()
            else:
                derep_assemblies = pd.concat([x[0] for x in dfs if None not in dfs])

            results_list = [x[1] for x in dfs if None not in dfs]
            if all(v is None for v in results_list):
                results = pd.DataFrame()
            else:
                results = pd.concat([x[1] for x in dfs if None not in dfs])

            stats_df_list = [x[3] for x in dfs if None not in dfs]
            if all(v is None for v in stats_df_list):
                stats_df = pd.DataFrame()
            else:
                stats_df = pd.concat([x[3] for x in dfs if None not in dfs])

            # Insert results into database
            insert_to_db(
                derep_assemblies=derep_assemblies,
                results=results,
                failed=failed,
                stats_df=stats_df,
                con=con,
                threads=args.threads,
                out_dir=args.out_dir,
                copy=args.copy,
            )

        jobs_done = retrieve_all_jobs_done(con)
        jobs_failed = retrieve_all_jobs_failed(con)
        genomes_derep = retrieve_all_genomes_derep(con)
        logging.info(
            "Jobs failed: {} / Jobs done: {}".format(
                jobs_failed.shape[0], jobs_done.shape[0]
            )
        )
        if not jobs_done.empty and not genomes_derep.empty:
            jobs_done = retrieve_all_jobs_done(con)
            genomes_derep = retrieve_all_genomes_derep(con)
            jobs_done = jobs_done.merge(genomes_derep)
            fname = prefix + "-derep-genomes_results.tsv"
            logging.info("Saving results to {}".format(fname))
            jobs_done.loc[:, "src"] = jobs_done["file"]
            if not args.copy:
                jobs_done.loc[:, "dst"] = ""
            else:
                jobs_done.loc[:, "dst"] = jobs_done["file"].apply(
                    lambda x: os.path.join(pathlib.Path(out_dir, pathlib.Path(x).name))
                )

                files = pd.DataFrame()
                files.loc[:, "src"] = jobs_done.loc[:, "src"]
                files.loc[:, "dst"] = jobs_done.loc[:, "dst"]
                files = list(files.itertuples(index=False, name=None))

                log.info(
                    "Copying dereplicated {:,} assemblies to output directory".format(
                        len(files)
                    )
                )

                if is_debug():
                    files = list(map(copy_files, files))
                else:
                    p = Pool(args.threads)
                    files = list(
                        tqdm(
                            p.imap_unordered(copy_files, files),
                            total=len(files),
                            leave=False,
                            ncols=80,
                        )
                    )
                p.close()
                p.join()
                log.info(
                    "Dereplication complete. Jobs saved to DB and files copied to {}\n".format(
                        out_dir
                    )
                )

            jobs_done["representative"] = jobs_done["representative"].astype(bool)
            jobs_done[["taxonomy", "accession", "representative", "src", "dst"]].to_csv(
                fname, index=False, sep="\t"
            )
    else:
        logging.info("Didn't find any new assemblies to process")

    con.close()


if __name__ == "__main__":
    main()
