import sqlite3
import pathlib
import urllib.parse
import logging
import os
from pathlib import Path
import pandas as pd

tables = [
    "taxa",
    "genomes",
    "genomes_derep",
    "results",
    "jobs_done",
    "jobs_failed",
    "stats",
]

log = logging.getLogger("my_logger")


def _path_to_uri(path):
    path = pathlib.Path(path)
    if path.is_absolute():
        return path.as_uri()
    return "file:" + urllib.parse.quote(path.as_posix(), safe=":/")


def check_if_db_exists(db):
    uri = _path_to_uri(db) + "?mode=rw"
    try:
        con = sqlite3.connect(uri, uri=True, isolation_level="EXCLUSIVE")
        log.info("Found DB in {}".format(db))

    except sqlite3.OperationalError:
        log.info("DB not found. Creating it")
        con = sqlite3.connect(db, isolation_level="EXCLUSIVE")
    return con


def check_if_db_empty(con):
    cursor = con.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    if cursor.fetchall():
        return False
    else:
        return True


def create_db_tables(con):
    """
    This function creates the tables needed to stores results and resuming failed jobs or updating new data
    """
    tax_table = "CREATE TABLE taxa(taxon TEXT PRIMARY KEY)"
    tax_table_idx = "CREATE INDEX idx_taxa_taxon ON taxa (taxon);"

    genomes_table = "CREATE TABLE genomes(taxon TEXT, accession TEXT PRIMARY KEY)"
    genomes_table_idx = "CREATE INDEX idx_genomes_taxon ON genomes (taxon);"

    genomes_derep_table = (
        "CREATE TABLE genomes_derep(accession TEXT PRIMARY KEY, representative INTEGER)"
    )
    genomes_derep_table_idx = (
        "CREATE INDEX idx_genomes_derep_acc ON genomes_derep (accession);"
    )

    results_table = "CREATE TABLE results(taxon TEXT PRIMARY KEY, weight REAL, communities INTEGER, n_genomes INTEGER, n_genomes_derep INTEGER)"
    results_table_idx = "CREATE INDEX idx_results_taxon ON results (taxon);"

    stats_table = "CREATE TABLE stats(taxon TEXT, representative TEXT PRIMARY KEY, n_nodes INTEGER, n_nodes_selected INTEGER, n_nodes_discarded INTEGER, graph_avg_weight REAL, graph_sd_weight REAL, graph_avg_weight_raw REAL, graph_sd_weight_raw REAL, subgraph_selected_avg_weight REAL, subgraph_selected_sd_weight REAL, subgraph_selected_avg_weight_raw REAL, subgraph_selected_sd_weight_raw REAL, subgraph_discarded_avg_weight REAL,subgraph_discarded_sd_weight REAL, subgraph_discarded_avg_weight_raw REAL, subgraph_discarded_sd_weight_raw REAL)"
    stats_table_idx = "CREATE INDEX idx_stats_taxon ON stats (representative);"

    jobs_done_table = (
        "CREATE TABLE jobs_done(taxon TEXT, accession TEXT PRIMARY KEY, file TEXT)"
    )
    jobs_done_table_idx = "CREATE INDEX idx_jobs_done_taxon ON jobs_done (taxon);"

    failed_table = "CREATE TABLE jobs_failed(taxon TEXT, accession TEXT PRIMARY KEY, file TEXT, reason TEXT)"
    failed_table_idx = "CREATE INDEX idx_jobs_failed_taxon ON jobs_failed (taxon);"

    # Create tables
    con.execute(tax_table)
    con.execute(genomes_table)
    con.execute(genomes_derep_table)
    con.execute(results_table)
    con.execute(jobs_done_table)
    con.execute(failed_table)
    con.execute(stats_table)
    # Create indices
    con.execute(tax_table_idx)
    con.execute(genomes_table_idx)
    con.execute(genomes_derep_table_idx)
    con.execute(results_table_idx)
    con.execute(failed_table_idx)
    con.execute(stats_table_idx)

    try:
        con.commit()
    except:
        pass


def check_db_tables(con):
    """
    This function check that the tables in the db are the correct ones
    """
    cursor = con.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    # tables = [
    #     "taxa",
    #     "genomes",
    #     "genomes_derep",
    #     "results",
    #     "jobs_done",
    # ]
    db_tables = cursor.fetchall()
    db_tables = [table[0] for table in db_tables]
    if set(tables) != set(db_tables):
        log.info("DB has incorrect tables. Dropping tables and recreating DB")
        for table in db_tables:
            con.execute("DROP TABLE %s" % table)
        try:
            con.commit()
        except:
            pass
        create_db_tables(con)
        log.info("DB has been recreated")
    else:
        log.info("DB tables seem to be OK")


def db_insert_taxa(con, taxon):
    """
    Function to insert in the DB the seen taxa
    """
    # "CREATE TABLE taxa(taxon TEXT PRIMARY KEY)"
    query = "INSERT INTO taxa (taxon) VALUES (?)"
    cursor = con.cursor()
    cursor.execute(query, (taxon,))


def db_insert_genomes(con, taxon, acc2assm):
    """
    Function to insert the genomes seen from each taxa and used for dereplication
    """
    # "CREATE TABLE genomes(taxon TEXT PRIMARY KEY, accession TEXT)"

    query = "INSERT INTO genomes (taxon, accession) VALUES (?, ?)"
    cursor = con.cursor()

    for assembly in acc2assm.keys():
        cursor.execute(query, (taxon, assembly))


def db_insert_genomes_derep(con, acc2assm, assms, reps):
    """
    Function to track the genomes that have been dereplicated
    """
    # "CREATE TABLE genomes_derep(accession TEXT PRIMARY KEY, representative INTEGER)"

    query = "INSERT INTO genomes_derep (accession, representative) VALUES (?, ?)"
    cursor = con.cursor()
    acc2assm = {k: v for k, v in acc2assm.items() if v in assms}

    for assembly in acc2assm.keys():
        if assembly in reps:
            is_rep = 1
        else:
            is_rep = 0
        cursor.execute(query, (assembly, is_rep))


def db_insert_results(con, taxon, weight, communities, n_genomes, n_genomes_derep):
    """
    Function to track the results from the dereplication
    """
    # "CREATE TABLE results(taxon TEXT PRIMARY KEY, weight REAL, communities INTEGER, n_genomes INTEGER, n_genomes_derep INTEGER)"
    query = "INSERT INTO results (taxon, weight, communities, n_genomes, n_genomes_derep) VALUES (?, ?, ?, ?, ?)"
    cursor = con.cursor()

    cursor.execute(query, (taxon, weight, communities, n_genomes, n_genomes_derep))


def db_insert_job_done(con, taxon, acc2assm, assms):
    """
    Function to keep track of the successful jobs and used for resuming a failed job
    """
    # "CREATE TABLE jobs_done(taxon TEXT, accession TEXT PRIMARY KEY, file TEXT)"
    query = "INSERT INTO jobs_done (taxon,accession,file) VALUES (?, ?, ?)"
    cursor = con.cursor()
    acc2assm = {k: v for k, v in acc2assm.items() if v in assms}

    for assembly in acc2assm.keys():
        cursor.execute(query, (taxon, assembly, os.path.basename(acc2assm[assembly])))


def check_if_done(con, taxon, acc2assm):
    # Find if taxon is in jobs done and taxa processed
    jobs_done = retrieve_jobs_done(con, taxon)
    genomes_done = retrieve_taxa_analyzed(con, taxon)

    if not jobs_done.empty and not genomes_done.empty:
        # check that there are no updates
        jobs_done = jobs_done.merge(genomes_done).loc[:, ["accession", "file"]]
        acc2assm.loc[:, "file"] = acc2assm["assembly"].apply(os.path.basename)

        needs_update = check_if_updates(
            jobs_done, acc2assm.loc[:, ["accession", "file"]]
        )
        if needs_update:
            remove_entries(taxon, tables, con)
            return False
        # Check files are in the folder
        # all_files_exist = not check_done_files_exists(files_done, out_dir)
        # if all_files_exist:
        #     return False
    else:
        return False
    return True


def remove_entries(taxon, tables, con):
    pass


def retrieve_jobs_done(con, taxon):
    query = "SELECT * from jobs_done where taxon =?"
    jobs = pd.read_sql(query, con, params=(taxon,))

    if not jobs.empty:
        return jobs
    else:
        return pd.DataFrame()


def retrieve_all_jobs_done(con):
    query = "SELECT * from jobs_done"
    jobs = pd.read_sql(query, con)

    if not jobs.empty:
        return jobs
    else:
        return pd.DataFrame()


def retrieve_all_jobs_failed(con):
    query = "SELECT * from jobs_failed"
    jobs = pd.read_sql(query, con)

    if not jobs.empty:
        return jobs
    else:
        return pd.DataFrame()


def retrieve_all_genomes_derep(con):
    query = "SELECT * from genomes_derep"
    jobs = pd.read_sql(query, con)

    if not jobs.empty:
        return jobs
    else:
        return pd.DataFrame()


def retrieve_taxa_analyzed(con, taxon):
    query = "SELECT * from genomes where taxon =?"
    taxons = pd.read_sql(query, con, params=(taxon,))

    if not taxons.empty:
        # accessions = [str(k[1]) for k in jobs]
        return taxons
    else:
        return pd.DataFrame()


def retrieve_all_taxa_analyzed(con):
    query = "SELECT * from taxa"
    taxons = pd.read_sql(query, con)

    if not taxons.empty:
        # accessions = [str(k[1]) for k in jobs]
        return taxons
    else:
        return pd.DataFrame()


# def retrieve_taxa_analyzed(con, taxon):
#     query = "SELECT * from genomes where taxon =?"
#     cursor = con.cursor()
#     cursor.execute(query, (taxon,))
#     jobs = cursor.fetchall()
#     if jobs:
#         accessions = [str(k[1]) for k in jobs]
#         return {taxon: accessions}
#     else:
#         return None


def retrieve_results_done(con, taxon):
    query = "SELECT * from results where taxon =?"
    cursor = con.cursor()
    cursor.execute(query, (taxon,))
    jobs = cursor.fetchall()
    if jobs:
        return True
    else:
        return False


def check_if_updates(old, new):
    df = pd.concat([new, old, old]).drop_duplicates(keep=False)
    if not df.empty:
        return True
    else:
        return False


def check_done_files_exists(files_done, out_dir):
    files = next(iter(files_done.values()))
    for file in files:
        file = Path(os.path.join(out_dir, file))
        if file.is_file():
            is_file = True
        else:
            return False
    return is_file


def delete_from_db(taxons, con):
    # first get accessions from the taxon of interest
    taxons = taxons["taxon"].tolist()

    placeholders = ", ".join(["?" for _ in taxons])
    query = "SELECT * FROM jobs_done WHERE taxon in ({});".format(placeholders)

    jobs_done = pd.read_sql_query(query, con, params=(taxons))

    cur = con.cursor()
    # Delete from taxa
    query = "DELETE FROM taxa WHERE taxon in ({});".format(placeholders)
    cur.execute(query, taxons)

    # Delete from jobs_done
    query = "DELETE FROM jobs_done WHERE taxon in ({});".format(placeholders)
    cur.execute(query, taxons)

    # Delete from genomes
    query = "DELETE FROM genomes WHERE taxon in ({});".format(placeholders)
    cur.execute(query, taxons)

    # Delete from results
    query = "DELETE FROM results WHERE taxon in ({});".format(placeholders)
    cur.execute(query, taxons)

    # Delete from genomes_derep
    accs = jobs_done["accession"].tolist()
    placeholders = ", ".join(["?" for _ in accs])
    query = "DELETE FROM genomes_derep WHERE accession in ({});".format(placeholders)

    cur.execute(query, accs)

    con.commit()
    return jobs_done
