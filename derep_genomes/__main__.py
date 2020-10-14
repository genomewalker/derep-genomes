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
    COPY,
    DEBUG,
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
)
from multiprocessing import Pool
import tqdm
from functools import partial
import pandas as pd


logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


def process_one_taxon(
    taxon,
    classification,
    out_dir,
    threads,
    threshold,
    chunks,
    slurm_config,
    tmp_dir,
    max_jobs_array,
    con,
):
    logging.info("Dereplicating {}".format(taxon))
    classification = classification[classification["taxon"] == taxon]

    acc_to_assemblies = classification.loc[:, ["accession", "assembly"]]
    n_assemblies = classification.shape[0]
    # check if already processed
    logging.info("Retrieving jobs done")
    is_done = check_if_done(con=con, taxon=taxon, acc2assm=acc_to_assemblies)

    if is_done:
        logging.info("Taxon already processed")
        return

    if n_assemblies == 0:
        return
    else:
        logging.info(
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

        logging.info("Saving data in DB")

        # Add taxa to taxa table
        derep_assemblies.loc[:, "taxon"] = taxon
        derep_assemblies.loc[:, "file"] = derep_assemblies["assembly"].apply(
            os.path.basename
        )

        results.loc[:, "taxon"] = taxon

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
            derep_assemblies.loc[:, ["accession", "representative"]],
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
        logging.info("Copying dereplicated assemblies to output directory")
        for assembly in derep_assemblies:
            if DEBUG:
                print("{} -> {}".format(assembly, out_dir))
            if COPY:
                files = pd.DataFrame()
                files.loc[:, "src"] = derep_assemblies.loc[:, "assembly"]
                files.loc[:, "dst"] = out_dir
                files = list(files.itertuples(index=False, name=None))
                if DEBUG:
                    files = list(map(copy_files, files))
                else:
                    p = Pool(threads)
                    files = list(
                        tqdm.tqdm(p.imap_unordered(copy_files, files), total=len(files))
                    )
        logging.info("Dereplication complete. Job saved in DB")


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
        logging.info("All entries in table {} already present".format(table))
    else:
        logging.info("Adding {} entries in table {}".format(query.shape[0], table))
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
        logging.info(
            "Copying {} files with singletons to {}".format(
                singletons.shape[0], out_dir
            )
        )
        files = pd.DataFrame()
        files.loc[:, "src"] = singletons.loc[:, "assembly"]
        files.loc[:, "dst"] = out_dir
        files = list(files.itertuples(index=False, name=None))

        if DEBUG:
            files = list(map(copy_files, files))
        else:
            p = Pool(threads)
            files = list(
                tqdm.tqdm(p.imap_unordered(copy_files, files), total=len(files))
            )


def find_assemblies(x, classifications, all_assm, debug):
    accessions = sorted(classifications[x])
    if debug:
        print(x)
        print(accessions[0])
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


def main():

    args = get_arguments()

    if COPY:
        out_dir = pathlib.Path(args.out_dir).absolute()
        os.makedirs(out_dir, exist_ok=True)

    tmp_dir = pathlib.Path(args.tmp_dir).absolute()

    con = check_if_db_exists(args.db)
    db_empty = check_if_db_empty(con)

    if db_empty:
        logging.info("Creating db tables")
        create_db_tables(con)
    else:
        logging.info("Checking correct tables exist")
        check_db_tables(con)

    all_assemblies = find_all_assemblies(args.in_dir)
    classifications = load_classifications(args.tax_file)

    logging.info("Filtering taxa with assemblies")
    classifications_df = pd.DataFrame(
        [
            (key, var, shorten_accession(var))
            for (key, L) in classifications.items()
            for var in L
        ],
        columns=["taxon", "accession", "accession_nover"],
    )

    all_assemblies_df = list(map(get_accession, all_assemblies))
    assm_data = classifications_df.merge(pd.DataFrame(all_assemblies_df))

    if args.selected_taxa:
        with args.selected_taxa as f:
            taxas = [line.rstrip() for line in f]
        assm_data_done = assm_data[assm_data["taxon"].isin(taxas)]
        assm_data = assm_data[~assm_data["taxon"].isin(taxas)]

    logging.info("Finding singletons...")
    taxon_counts = assm_data.groupby("taxon", as_index=False)["taxon"].agg(
        {"count": "count"}
    )

    classification_singletons = taxon_counts[taxon_counts["count"] < 2][
        "taxon"
    ].tolist()
    classification_singletons = assm_data[
        assm_data["taxon"].isin(classification_singletons)
    ].reset_index(drop=True)

    classifications_sorted_counts = taxon_counts[taxon_counts["count"] > 1].sort_values(
        by=["count"], ascending=True
    )
    classifications_sorted = assm_data[
        assm_data["taxon"].isin(classifications_sorted_counts["taxon"])
    ].reset_index(drop=True)
    logging.info("Processing singletons")
    process_sigletons(
        singletons=classification_singletons,
        out_dir=args.out_dir,
        threads=args.threads,
        con=con,
    )

    # derep_assemblies, results, reps = dereplicate(
    #     acc_to_assemblies,
    #     threads,
    #     threshold,
    #     chunks,
    #     slurm_config,
    #     tmp_dir,
    #     debug,
    #     max_jobs_array,
    #     con,
    # )
    taxons = classifications_sorted_counts.loc[:, ["taxon"]]

    for taxa in taxons["taxon"]:
        results = process_one_taxon(
            taxon=taxa,
            classification=classifications_sorted,
            out_dir=args.out_dir,
            threads=args.threads,
            threshold=args.threshold,
            chunks=args.chunks,
            slurm_config=args.slurm_config,
            tmp_dir=tmp_dir,
            max_jobs_array=args.max_jobs_array,
            con=con,
        )

    # for taxon in classifications_sorted:
    #     process_one_taxon(
    #         classification=taxon,
    #         accessions=classifications[taxon],
    #         out_dir=out_dir,
    #         threads=args.threads,
    #         threshold=args.threshold,
    #         chunks=args.chunks,
    #         slurm_config=args.slurm_config,
    #         tmp_dir=tmp_dir,
    #         debug=args.debug,
    #         max_jobs_array=args.max_jobs_array,
    #         con=con,
    #     )
    con.close()


if __name__ == "__main__":
    main()
