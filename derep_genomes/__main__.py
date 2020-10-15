"""
This script will dereplicate GTDB assemblies to a user-specified threshold, copying the
dereplicated assemblies to a new directory.

Usage:
    dereplicate_assemblies.py --threshold 0.005 assemblies derep bac_and_arc_taxonomy_r86.tsv

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

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
)
from derep_genomes.graph import dereplicate
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
)
from multiprocessing import Pool
import tqdm
from functools import partial
import pandas as pd

COPY = False

log = logging.getLogger("my_logger")


def process_one_taxon(taxon, parms):
    classification = parms["classification"]
    threads = parms["threads"]
    threshold = parms["threshold"]
    chunks = parms["chunks"]
    slurm_config = parms["slurm_config"]
    tmp_dir = parms["tmp_dir"]
    max_jobs_array = parms["max_jobs_array"]

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
        derep_assemblies, results = dereplicate(
            all_assemblies=acc_to_assemblies,
            threads=threads,
            threshold=threshold,
            chunks=chunks,
            slurm_config=slurm_config,
            tmp_dir=tmp_dir,
            max_jobs_array=max_jobs_array,
        )

        # Add taxa to taxa table
        derep_assemblies.loc[:, "taxon"] = taxon
        derep_assemblies.loc[:, "file"] = derep_assemblies["assembly"].apply(
            os.path.basename
        )
        results.loc[:, "taxon"] = taxon
        return derep_assemblies, results


def insert_to_db(derep_assemblies, results, con, threads, out_dir):
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

    # Add jobs done to table
    query = check_existing(
        derep_assemblies.loc[:, ["taxon", "accession", "file"]],
        columns=["taxon", "accession", "file"],
        table="jobs_done",
        con=con,
    )
    jobs_done_db = insert_pd_sql(query=query, table="jobs_done", con=con)

    try:
        con.commit()
    except:
        pass
    for assembly in derep_assemblies:
        log.debug("{} -> {}".format(assembly, out_dir))
        if COPY:
            log.info("Copying dereplicated assemblies to output directory")
            files = pd.DataFrame()
            files.loc[:, "src"] = derep_assemblies.loc[:, "assembly"]
            files.loc[:, "dst"] = out_dir
            files = list(files.itertuples(index=False, name=None))
            if is_debug():
                files = list(map(copy_files, files))
            else:
                p = Pool(threads)
                files = list(
                    tqdm.tqdm(p.imap_unordered(copy_files, files), total=len(files))
                )
    log.info("Dereplication complete. Jobs saved to DB\n")


def copy_files(src_dst):
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


def process_sigletons(singletons, out_dir, threads, con):

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

    if COPY:
        log.info(
            "Copying {} files with singletons to {}".format(
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
                tqdm.tqdm(p.imap_unordered(copy_files, files), total=len(files))
            )


def find_assemblies(x, classifications, all_assm):
    accessions = sorted(classifications[x])
    log.debug(x)
    log.debug(accessions[0])
    res = find_assemblies_for_accessions(accessions=accessions, all_assemblies=all_assm)
    return res


def shorten_accession(accession):
    if accession.startswith("GCF_") or accession.startswith("GCA_"):
        accession = accession.split(".")[0]
        assert len(accession) == 13
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

    COPY = args.copy

    if COPY:
        out_dir = pathlib.Path(args.out_dir).absolute()
        os.makedirs(out_dir, exist_ok=True)

    tmp_dir = pathlib.Path(args.tmp_dir).absolute()

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
    log.info("Found {:,} assemblies".format(classifications_df.shape[0]))
    all_assemblies_df = list(map(get_accession, all_assemblies))
    assm_data = classifications_df.merge(pd.DataFrame(all_assemblies_df))
    assm_data.loc[:, "file"] = assm_data["assembly"].apply(os.path.basename)

    if args.selected_taxa:
        with args.selected_taxa as f:
            taxas = [line.rstrip() for line in f]
        assm_data = assm_data[assm_data["taxon"].isin(taxas)]

    taxa_in_db = retrieve_all_taxa_analyzed(con)
    jobs_done = retrieve_all_jobs_done(con)
    ntaxa = len(list(set(assm_data["taxon"])))

    log.info("Found {:,}/{:,} taxa in the DB".format(taxa_in_db.shape[0], ntaxa))

    log.info("Checking assemblies already processed")

    if not jobs_done.empty and not taxa_in_db.empty:
        # check that there are no updates
        jobs_done = jobs_done.merge(taxa_in_db).loc[:, ["accession", "file"]]
        to_do = pd.concat(
            [assm_data.loc[:, ["accession", "file"]], jobs_done, jobs_done]
        ).drop_duplicates(keep=False)
        to_do = assm_data[assm_data["accession"].isin(to_do["accession"])]
        to_rm = pd.concat(
            [to_do.loc[:, ["taxon"]], taxa_in_db, taxa_in_db]
        ).drop_duplicates(keep=False)
        # TODO: create funtion that removes entries that need to be updated
        log.info(
            "{:,}/{:,} assemblies already processed".format(
                to_do.shape[0], assm_data.shape[0]
            )
        )
    else:
        log.info("No jobs found in the database")
        to_do = assm_data

    log.info("Finding singletons...")
    taxon_counts = to_do.groupby("taxon", as_index=False)["taxon"].agg(
        {"count": "count"}
    )

    classification_singletons = taxon_counts[taxon_counts["count"] < 2][
        "taxon"
    ].tolist()
    classification_singletons = to_do[
        to_do["taxon"].isin(classification_singletons)
    ].reset_index(drop=True)

    # classifications_sorted_counts = taxon_counts[taxon_counts["count"] > 1].sort_values(
    #     by=["count"], ascending=True
    # )

    assm_min = 2
    assm_max = 5

    classification_small = taxon_counts[
        taxon_counts["count"].isin(range(assm_min, assm_max + 1))
    ]["taxon"].tolist()
    classification_small = to_do[to_do["taxon"].isin(classification_small)].reset_index(
        drop=True
    )

    classification_large = taxon_counts[taxon_counts["count"] > 5]["taxon"].tolist()
    classification_large = to_do[to_do["taxon"].isin(classification_large)].reset_index(
        drop=True
    )

    log.info("Processing {:,} singletons".format(classification_singletons.shape[0]))
    if not classification_singletons.empty:
        process_sigletons(
            singletons=classification_singletons,
            out_dir=args.out_dir,
            threads=args.threads,
            con=con,
        )
    else:
        log.info("No singletons found\n")

    if not classification_small.empty:
        taxons = list(set(classification_small["taxon"]))

        parms = {
            "classification": classification_small,
            "threads": 2,
            "threshold": args.threshold,
            "chunks": args.chunks,
            "slurm_config": None,
            "tmp_dir": tmp_dir,
            "max_jobs_array": args.max_jobs_array
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
            files = list(map(partial(process_one_taxon, parms=parms), taxons))
        else:
            p = Pool(nworkers)
            dfs = list(
                tqdm.tqdm(
                    p.imap_unordered(partial(process_one_taxon, parms=parms), taxons),
                    total=len(taxons),
                )
            )

        derep_assemblies = pd.concat([x[0] for x in dfs])
        results = pd.concat([x[1] for x in dfs])
        insert_to_db(
            derep_assemblies=derep_assemblies,
            results=results,
            con=con,
            threads=args.threads,
            out_dir=args.out_dir,
        )

    if not classification_large.empty:
        taxons = list(set(classification_large["taxon"]))
        log.info(
            "Dereplicating {:,} taxa with more assemblies using SLURM".format(
                len(taxons), assm_min, assm_max
            )
        )
        parms_large = {
            "classification": classification_large,
            "threads": args.threads,
            "threshold": args.threshold,
            "chunks": args.chunks,
            "slurm_config": args.slurm_config.name,
            "tmp_dir": tmp_dir,
            "max_jobs_array": args.max_jobs_array
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
                )
            )
        derep_assemblies = pd.concat([x[0] for x in dfs])
        results = pd.concat([x[1] for x in dfs])
        insert_to_db(
            derep_assemblies=derep_assemblies,
            results=results,
            con=con,
            threads=args.threads,
            out_dir=args.out_dir,
        )

    con.close()


if __name__ == "__main__":
    main()
