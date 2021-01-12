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

from scipy import stats
from derep_genomes import __version__

import subprocess
import sys
import os
import shutil
from derep_genomes.general import (
    get_arguments,
    find_all_assemblies,
    load_classifications,
    find_assemblies_for_accessions,
    suppress_stdout,
    get_assembly_filename,
    applyParallel,
    is_debug,
    get_assembly_length,
)
from derep_genomes.graph import dereplicate_ANI, dereplicate_mash, concat_df
import logging
import pathlib
import sqlite3
from derep_genomes.dbops import (
    check_if_db_exists,
    check_if_db_empty,
    create_db_tables,
    check_db_tables,
    retrieve_jobs_done,
    db_insert_job_done,
    db_insert_taxa,
    db_insert_genomes,
    db_insert_genomes_derep,
    db_insert_results,
    check_if_done,
    retrieve_all_taxa_analyzed,
    retrieve_all_jobs_done,
    retrieve_all_jobs_failed,
    delete_from_db,
    retrieve_all_genomes_derep,
)
from multiprocessing import Pool
import tqdm
from functools import partial
import pandas as pd
import time
import tempfile

log = logging.getLogger("my_logger")
timestr = time.strftime("%Y%m%d-%H%M%S")


def process_one_taxon(taxon, parms):
    classification = parms["classification"]
    threads = parms["threads"]
    threshold = parms["threshold"]
    chunks = parms["chunks"]
    slurm_config = parms["slurm_config"]
    tmp_dir = parms["tmp_dir"]
    max_jobs_array = parms["max_jobs_array"]
    mash_threshold = parms["mash_threshold"]
    failed_df = None

    log.debug("Dereplicating {}".format(taxon))
    classification = classification[classification["taxon"] == taxon]
    acc_to_assemblies = classification.loc[:, ["accession", "assembly"]]
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

        log.debug("Preclustering with MASH...")
        acc_to_assemblies_mash = dereplicate_mash(
            all_assemblies=acc_to_assemblies,
            threads=threads,
            tmp_dir=tmp_dir,
            mash_threshold=mash_threshold,
            threshold=threshold,
        )
        acc_to_assemblies_mash_r = acc_to_assemblies_mash[
            acc_to_assemblies_mash["representative"] == 1
        ].copy()
        acc_to_assemblies_mash_d = acc_to_assemblies_mash[
            acc_to_assemblies_mash["derep"] == 1
        ].copy()
        acc_to_assemblies_mash_o = acc_to_assemblies_mash[
            acc_to_assemblies_mash["derep"] == 0
        ].copy()

        log.debug(
            "Identified {:,} linked assemblies in {:,} MASH communities".format(
                acc_to_assemblies_mash_d.shape[0], acc_to_assemblies_mash_r.shape[0]
            )
        )

        # We refine each MASH community
        log.debug(
            "Refining {:,} MASH selected assemblies".format(
                acc_to_assemblies_mash_d.shape[0]
            )
        )
        derep_assemblies, results, failed, stats_df = dereplicate_ANI(
            all_assemblies=acc_to_assemblies_mash_d,
            threads=threads,
            threshold=threshold,
            chunks=chunks,
            slurm_config=slurm_config,
            tmp_dir=tmp_dir,
            max_jobs_array=max_jobs_array,
        )
        if failed:
            failed_df = acc_to_assemblies.copy()
            failed_df.loc[:, "file"] = failed_df["assembly"].apply(os.path.basename)
            failed_df.loc[:, "taxon"] = taxon
            failed_df.loc[:, "reason"] = failed
            failed_df = failed_df[["taxon", "accession", "file", "reason"]]
        else:
            # Add taxa to taxa table
            derep_assemblies = concat_df([derep_assemblies, acc_to_assemblies_mash_o])
            derep_assemblies.loc[:, "taxon"] = taxon
            derep_assemblies.loc[:, "file"] = derep_assemblies["assembly"].apply(
                os.path.basename
            )
            results.loc[:, "taxon"] = taxon
            if stats_df is not None:
                stats_df.loc[:, "taxon"] = taxon

        return derep_assemblies, results, failed_df, stats_df


def insert_to_db(
    derep_assemblies, results, failed, stats_df, con, threads, out_dir, copy
):
    if not derep_assemblies.empty:
        query = check_existing(
            df=derep_assemblies.loc[:, ["taxon"]].drop_duplicates(),
            columns=["taxon"],
            table="taxa",
            con=con,
        )
        insert_pd_sql(query=query, table="taxa", con=con)

        # Add genomes to genome table
        query = check_existing(
            df=derep_assemblies.loc[:, ["taxon", "accession"]],
            columns=["taxon", "accession"],
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
                :, ["taxon", "weight", "communities", "n_genomes", "n_genomes_derep"]
            ],
            columns=["taxon", "weight", "communities", "n_genomes", "n_genomes_derep"],
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
                    "taxon",
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
                "taxon",
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
            derep_assemblies.loc[:, ["taxon", "accession", "file"]],
            columns=["taxon", "accession", "file"],
            table="jobs_done",
            con=con,
        )
        jobs_done_db = insert_pd_sql(query=query, table="jobs_done", con=con)

    # Add failed jobs to table
    if not failed.empty:
        query = check_existing(
            failed.loc[:, ["taxon", "accession", "file", "reason"]],
            columns=["taxon", "accession", "file", "reason"],
            table="jobs_failed",
            con=con,
        )
        jobs_failed = insert_pd_sql(query=query, table="jobs_failed", con=con)

        query = check_existing(
            df=failed.loc[:, ["taxon"]].drop_duplicates(),
            columns=["taxon"],
            table="taxa",
            con=con,
        )
        insert_pd_sql(query=query, table="taxa", con=con)

        query = check_existing(
            failed.loc[:, ["taxon", "accession", "file"]],
            columns=["taxon", "accession", "file"],
            table="jobs_done",
            con=con,
        )
        jobs_failed = insert_pd_sql(query=query, table="jobs_done", con=con)

    try:
        con.commit()
    except:
        pass
    if copy:
        files = pd.DataFrame()
        files.loc[:, "src"] = derep_assemblies[derep_assemblies["derep"] == 1].loc[
            :, "assembly"
        ]
        files.loc[:, "dst"] = out_dir
        files = list(files.itertuples(index=False, name=None))
        log.info(
            "Copying dereplicated {:,} assemblies to output directory".format(
                len(files)
            )
        )

        if is_debug():
            files = list(map(copy_files, files))
        else:
            p = Pool(threads)
            files = list(
                tqdm.tqdm(
                    p.imap_unordered(copy_files, files),
                    total=len(files),
                    leave=False,
                    ncols=80,
                )
            )

        log.info(
            "Dereplication complete. Jobs saved to DB and files copied to {}\n".format(
                out_dir
            )
        )
        return
    log.info("Dereplication complete. Jobs saved to DB\n")


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

    df = pd.concat([new, old, old]).drop_duplicates(keep=False)
    return df


def insert_pd_sql(query, table, con):
    if query.empty:
        log.debug("All entries in table {} already present".format(table))
    else:
        log.debug("Adding {} entries in table {}".format(query.shape[0], table))
        query.to_sql(name=table, con=con, if_exists="append", index=False)


def process_sigletons(singletons, out_dir, threads, con, copy):

    # Insert data to table taxa
    query = check_existing(singletons, columns=["taxon"], table="taxa", con=con)
    insert_pd_sql(query=query, table="taxa", con=con)

    # Insert data to table genomes
    query = check_existing(
        singletons, columns=["taxon", "accession"], table="genomes", con=con
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
        columns=["taxon", "weight", "communities", "n_genomes", "n_genomes_derep"],
        table="results",
        con=con,
    )
    results = insert_pd_sql(query=query, table="results", con=con)
    # Insert data to table jobs_done
    singletons.loc[:, "file"] = singletons["assembly"].apply(
        lambda x: os.path.basename(x)
    )
    query = check_existing(
        singletons, columns=["taxon", "accession", "file"], table="jobs_done", con=con
    )
    jobs_done = insert_pd_sql(query=query, table="jobs_done", con=con)

    try:
        con.commit()
    except:
        pass

    if copy:
        log.info(
            "Copying {} files with singletons to {}\n".format(
                singletons.shape[0], out_dir
            )
        )
        files = pd.DataFrame()
        files.loc[:, "src"] = singletons.loc[:, "assembly"]
        files.loc[:, "dst"] = out_dir
        files = list(files.itertuples(index=False, name=None))

        if is_debug():
            files = list(map(copy_files, files))
        else:
            p = Pool(threads)
            files = list(
                tqdm.tqdm(
                    p.imap_unordered(copy_files, files),
                    total=len(files),
                    leave=False,
                    ncols=80,
                )
            )


def find_assemblies(x, classifications, all_assm):
    accessions = sorted(classifications[x])
    log.debug(x)
    log.debug(accessions[0])
    res = find_assemblies_for_accessions(accessions=accessions, all_assemblies=all_assm)
    return res

#TODO: Improve the matching between filenames and accessions
def shorten_accession(accession):
    if accession.startswith("GCF_") or accession.startswith("GCA_"):
        accession = accession.split(".")[0]
        assert len(accession) == 13
    else:
        accession = accession.split(".")[0]
    return accession


def get_accession(assembly):
    res = {}
    accession = os.path.basename(assembly)
    accession = shorten_accession(accession)
    res = {"accession_nover": accession, "assembly": assembly}
    return res


def check_done(data, con):
    data.loc[:, "file"] = data["assembly"].apply(os.path.basename)
    df = check_existing(
        data, columns=["taxon", "accession", "file"], table="jobs_done", con=con
    )
    return df


def check_done_apply(df, parms):
    db = parms["db"]
    con = sqlite3.connect(db, isolation_level="EXCLUSIVE")
    taxon = list(set(df["taxon"]))[0]
    acc_to_assemblies = df.loc[:, ["accession", "assembly"]]
    is_done = check_if_done(con=con, taxon=taxon, acc2assm=acc_to_assemblies)

    con.close()
    return is_done


def mute():
    sys.stdout = open(os.devnull, "w")


def main():

    logging.basicConfig(
        level=logging.DEBUG, format="%(levelname)s ::: %(asctime)s ::: %(message)s"
    )

    args = get_arguments()
    logging.getLogger("my_logger").setLevel(
        logging.DEBUG if args.debug else logging.INFO
    )

    if args.copy:
        out_dir = pathlib.Path(args.out_dir).absolute()
        os.makedirs(out_dir, exist_ok=True)

    in_dir = pathlib.Path(args.in_dir).absolute()

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

    all_assemblies = find_all_assemblies(args.in_dir)
    classifications = load_classifications(args.tax_file)

    log.info("Filtering taxa with assemblies")
    classifications_df = pd.DataFrame(
        [
            (key, var, shorten_accession(var))
            for (key, L) in classifications.items()
            for var in L
        ],
        columns=["taxon", "accession", "accession_nover"],
    )
    log.info(
        "Found {:,} assembly files in {}".format(
            classifications_df.shape[0], args.in_dir
        )
    )
    all_assemblies_df = list(map(get_accession, all_assemblies))
    assm_data = classifications_df.merge(pd.DataFrame(all_assemblies_df))

    assm_data.loc[:, "file"] = assm_data["assembly"].apply(os.path.basename)

    if args.selected_taxa:
        with args.selected_taxa as f:
            taxas = [line.rstrip() for line in f]
        assm_data = assm_data[assm_data["taxon"].isin(taxas)]

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
        to_rm = pd.concat([to_do.loc[:, ["taxon"]].drop_duplicates(), taxa_in_db])
        to_rm = to_rm[to_rm["taxon"].duplicated()].reset_index()

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

            rfname = timestr + "-derep-genomes_removed-db-entries.tsv"
            removed.to_csv(rfname, index=False, sep="\t")

    else:
        log.info("No jobs found in the database\n")
        to_do = assm_data

    if not to_do.empty:
        assm_min = 2
        assm_max = 10

        # log.info("Calculating assembly lengths")

        # p = Pool(args.threads)
        # to_do.loc[:, "len"] = list(
        #     tqdm.tqdm(
        #         p.imap(get_assembly_length, to_do["assembly"].to_list()),
        #         total=to_do.shape[0],
        #         leave=False,
        #         ncols=80,
        #     )
        # )
        # print(to_do)
        # exit(0)

        log.info(
            "Splitting taxa in groups of singleton, between {}-{} and more than {} assemblies\n".format(
                assm_min, assm_max, assm_max
            )
        )

        taxon_counts = to_do.groupby("taxon", as_index=False)["taxon"].agg(
            {"count": "count"}
        )

        classification_singletons = taxon_counts[taxon_counts["count"] < 2][
            "taxon"
        ].tolist()
        classification_singletons = to_do[
            to_do["taxon"].isin(classification_singletons)
        ].reset_index(drop=True)

        classification_small = taxon_counts[
            taxon_counts["count"].isin(range(assm_min, assm_max + 1))
        ]["taxon"].tolist()
        classification_small = to_do[
            to_do["taxon"].isin(classification_small)
        ].reset_index(drop=True)

        classification_large = taxon_counts[taxon_counts["count"] > assm_max][
            "taxon"
        ].tolist()
        classification_large = to_do[
            to_do["taxon"].isin(classification_large)
        ].reset_index(drop=True)

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
        # else:
        #     log.info("No singletons found\n")

        if not classification_small.empty:
            taxons = list(set(classification_small["taxon"]))

            parms = {
                "classification": classification_small,
                "threads": 2,
                "threshold": args.threshold,
                "chunks": args.chunks,
                "slurm_config": None,
                "tmp_dir": tmp_dir,
                "max_jobs_array": args.max_jobs_array,
                "mash_threshold": args.mash_threshold,
            }

            if args.threads > 2:
                nworkers = round(args.threads / 2)
            else:
                nworkers = 1
            log.info(
                "Dereplicating {:,} taxa with {} to {} assemblies using {} workers".format(
                    len(taxons), assm_min, assm_max, nworkers
                )
            )

            if is_debug():
                dfs = list(map(partial(process_one_taxon, parms=parms), taxons))
            else:
                p = Pool(nworkers)
                dfs = list(
                    tqdm.tqdm(
                        p.imap_unordered(
                            partial(process_one_taxon, parms=parms), taxons
                        ),
                        total=len(taxons),
                        leave=False,
                        ncols=80,
                    )
                )

            l = [x[2] for x in dfs if None not in dfs]
            if all(v is None for v in l):
                failed = pd.DataFrame()
            else:
                failed = pd.concat([x[2] for x in dfs if not None in dfs])

            l = [x[0] for x in dfs if None not in dfs]
            if all(v is None for v in l):
                derep_assemblies = pd.DataFrame()
            else:
                derep_assemblies = pd.concat([x[0] for x in dfs if not None in dfs])

            l = [x[1] for x in dfs if None not in dfs]
            if all(v is None for v in l):
                results = pd.DataFrame()
            else:
                results = pd.concat([x[1] for x in dfs if not None in dfs])
            if all(v is None for v in l):
                stats_df = pd.DataFrame()
            else:
                stats_df = pd.concat([x[3] for x in dfs if not None in dfs])

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

        if not classification_large.empty:
            taxons = list(set(classification_large["taxon"]))
            if args.slurm_config is not None:
                log.info(
                    "Dereplicating {:,} taxa with more than {} assemblies using SLURM".format(
                        len(taxons), assm_max
                    )
                )
            else:
                log.info(
                    "Dereplicating {:,} taxa with more than 5 assemblies using {} threads".format(
                        len(taxons), args.threads
                    )
                )
                log.warning("This can take a long time!!!")
            parms_large = {
                "classification": classification_large,
                "threads": args.threads,
                "threshold": args.threshold,
                "chunks": args.chunks,
                "slurm_config": args.slurm_config.name,
                "tmp_dir": tmp_dir,
                "max_jobs_array": args.max_jobs_array,
                "mash_threshold": args.mash_threshold,
            }

            if is_debug():
                dfs = list(
                    map(partial(process_one_taxon, parms=parms_large), taxons),
                )
            else:
                dfs = list(
                    tqdm.tqdm(
                        map(partial(process_one_taxon, parms=parms_large), taxons),
                        total=len(taxons),
                        leave=False,
                        ncols=80,
                    )
                )

            l = [x[2] for x in dfs if None not in dfs]
            if all(v is None for v in l):
                failed = pd.DataFrame()
            else:
                failed = pd.concat([x[2] for x in dfs if not None in dfs])

            l = [x[0] for x in dfs if None not in dfs]
            if all(v is None for v in l):
                derep_assemblies = pd.DataFrame()
            else:
                derep_assemblies = pd.concat([x[0] for x in dfs if not None in dfs])

            l = [x[1] for x in dfs if None not in dfs]
            if all(v is None for v in l):
                results = pd.DataFrame()
            else:
                results = pd.concat([x[1] for x in dfs if not None in dfs])
            if all(v is None for v in l):
                stats_df = pd.DataFrame()
            else:
                stats_df = pd.concat([x[3] for x in dfs if not None in dfs])
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
            fname = timestr + "-derep-genomes_results.tsv"
            logging.info("Saving results to {}".format(fname))
            jobs_done.loc[:, "src"] = jobs_done["file"].apply(
                lambda x: os.path.join(in_dir, x)
            )
            if not args.copy:
                jobs_done.loc[:, "dst"] = ""
            else:
                jobs_done.loc[:, "dst"] = jobs_done["file"].apply(
                    lambda x: os.path.join(out_dir, x)
                )
            jobs_done["representative"] = jobs_done["representative"].astype(bool)
            jobs_done[["taxon", "accession", "representative", "src", "dst"]].to_csv(
                fname, index=False, sep="\t"
            )
    else:
        logging.info("Didn't find any new assemblies to process")
    con.close()


if __name__ == "__main__":
    main()
