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

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


def process_one_taxon(
    classification,
    accessions,
    all_assemblies,
    out_dir,
    threads,
    threshold,
    chunks,
    slurm_config,
    tmp_dir,
    debug,
    max_jobs_array,
    con,
):
    accessions = sorted(accessions)
    print()
    logging.info("Dereplicating {}".format(classification))
    acc_to_assemblies = find_assemblies_for_accessions(accessions, all_assemblies)

    # check if already processed
    logging.info("Retrieving jobs done")
    is_done = check_if_done(
        con=con, taxon=classification, acc2assm=acc_to_assemblies, out_dir=out_dir
    )

    if is_done:
        logging.info("Taxon already processed")
        return

    if len(acc_to_assemblies) == 0:
        return
    if len(acc_to_assemblies) == 1:
        rep = str(next(iter(acc_to_assemblies.keys())))
        only_assembly = list(acc_to_assemblies.values())[0]

        logging.info("Only one assembly for this species, copying to output directory:")
        if debug:
            print("{} -> {}".format(only_assembly, out_dir))
        shutil.copy(only_assembly, out_dir)
        # add taxa
        # add accessions
        # add results
        # add job done

        logging.info("Saving data in DB")
        db_insert_taxa(con=con, taxon=classification)
        db_insert_genomes(con=con, taxon=classification, acc2assm=acc_to_assemblies)
        db_insert_genomes_derep(
            con=con, acc2assm=acc_to_assemblies, assms=only_assembly, reps=rep
        )
        db_insert_results(
            con=con,
            taxon=classification,
            weight=0,
            communities=1,
            n_genomes=1,
            n_genomes_derep=1,
        )
        db_insert_job_done(
            con=con,
            taxon=classification,
            acc2assm=acc_to_assemblies,
            assms=only_assembly,
        )
        try:
            con.commit()
        except:
            pass
        logging.info("Dereplication complete. Job saved in DB")

    else:
        logging.info(
            "{:,} assemblies for this species, clustering to dereplicate.".format(
                len(acc_to_assemblies)
            )
        )
        derep_assemblies, results, reps = dereplicate(
            acc_to_assemblies,
            threads,
            threshold,
            chunks,
            slurm_config,
            tmp_dir,
            debug,
            max_jobs_array,
            con,
        )
        logging.info("Copying dereplicated assemblies to output directory")
        for assembly in derep_assemblies:
            if debug:
                print("{} -> {}".format(assembly, out_dir))
            shutil.copy(assembly, out_dir)

        logging.info("Saving data in DB")
        db_insert_taxa(con=con, taxon=classification)
        db_insert_genomes(con=con, taxon=classification, acc2assm=acc_to_assemblies)
        db_insert_genomes_derep(
            con=con, acc2assm=acc_to_assemblies, assms=derep_assemblies, reps=reps
        )
        db_insert_results(
            con=con,
            taxon=classification,
            weight=results[0],
            communities=results[1],
            n_genomes=results[2],
            n_genomes_derep=results[3],
        )
        db_insert_job_done(
            con=con,
            taxon=classification,
            acc2assm=acc_to_assemblies,
            assms=derep_assemblies,
        )
        try:
            con.commit()
        except:
            pass
        logging.info("Dereplication complete. Job saved in DB")


def main():
    args = get_arguments()
    out_dir = pathlib.Path(args.out_dir).absolute()
    os.makedirs(out_dir, exist_ok=True)

    db = os.path.join(out_dir, "derep-genomes.db")

    con = check_if_db_exists(db)
    db_empty = check_if_db_empty(con)

    if db_empty:
        logging.info("Creating db tables")
        create_db_tables(con)
    else:
        logging.info("Checking correct tables exist")
        check_db_tables(con)

    all_assemblies = find_all_assemblies(args.in_dir)
    classifications = load_classifications(args.tax_file)

    if args.selected_taxa:
        with args.selected_taxa as f:
            taxas = [line.rstrip() for line in f]
        classifications = {k: classifications[k] for k in taxas if k in classifications}

    classifications_sorted = [
        k
        for k, v in sorted(
            classifications.items(), key=lambda item: len(item[1]), reverse=False
        )
    ]

    tmp_dir = pathlib.Path(args.tmp_dir).absolute()

    for taxon in classifications_sorted:
        process_one_taxon(
            classification=taxon,
            accessions=classifications[taxon],
            all_assemblies=all_assemblies,
            out_dir=out_dir,
            threads=args.threads,
            threshold=args.threshold,
            chunks=args.chunks,
            slurm_config=args.slurm_config,
            tmp_dir=tmp_dir,
            debug=args.debug,
            max_jobs_array=args.max_jobs_array,
            con=con,
        )
    con.close()


if __name__ == "__main__":
    main()
