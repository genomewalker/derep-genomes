import sqlite3
import pathlib
import urllib.parse
import logging
import os
from pathlib import Path

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


def _path_to_uri(path):
    path = pathlib.Path(path)
    if path.is_absolute():
        return path.as_uri()
    return "file:" + urllib.parse.quote(path.as_posix(), safe=":/")


def check_if_db_exists(db):
    uri = _path_to_uri(db) + "?mode=rw"
    try:
        con = sqlite3.connect(uri, uri=True, isolation_level="EXCLUSIVE")
        logging.info("Found DB in {}".format(db))

    except sqlite3.OperationalError:
        logging.info("DB not found. Creating it")
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
    genomes_table = "CREATE TABLE genomes(taxon TEXT, accession TEXT PRIMARY KEY)"
    derep_genomes_table = (
        "CREATE TABLE genomes_derep(accession TEXT PRIMARY KEY, representative INTEGER)"
    )
    results_table = "CREATE TABLE results(taxon TEXT PRIMARY KEY, weight REAL, communities INTEGER, n_genomes INTEGER, n_genomes_derep INTEGER)"
    jobs_done_table = (
        "CREATE TABLE jobs_done(taxon TEXT, accession TEXT PRIMARY KEY, file TEXT)"
    )
    con.execute(tax_table)
    con.execute(genomes_table)
    con.execute(derep_genomes_table)
    con.execute(results_table)
    con.execute(jobs_done_table)
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
    tables = [
        "taxa",
        "genomes",
        "genomes_derep",
        "results",
        "jobs_done",
    ]
    db_tables = cursor.fetchall()
    db_tables = [table[0] for table in db_tables]
    if set(tables) != set(db_tables):
        logging.info("DB has incorrect tables. Dropping tables and recreating DB")
        for table in db_tables:
            con.execute("DROP TABLE %s" % table)
        try:
            con.commit()
        except:
            pass
        create_db_tables(con)
        logging.info("DB has been recreated")
    else:
        logging.info("DB tables seem to be OK")


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


def check_if_done(con, taxon, acc2assm, out_dir):
    # Find if taxon is in jobs done and taxa processed
    jobs_done, files_done = retrieve_jobs_done(con, taxon)
    genomes_done = retrieve_taxa_analyzed(con, taxon)
    # genomes_derep_done = retrieve_taxa_derep(con, taxon)
    # results_done = retrieve_results_done(con, taxon)
    if jobs_done and genomes_done:
        # check that there are no updates
        needs_update = check_if_updates(genomes_done, acc2assm)
        if needs_update:
            return False
        # Check files are in the folder
        all_files_exist = not check_done_files_exists(files_done, out_dir)
        if all_files_exist:
            return False
    else:
        return False

    return True


def retrieve_jobs_done(con, taxon):
    query = "SELECT * from jobs_done where taxon =?"
    cursor = con.cursor()
    cursor.execute(query, (taxon,))
    jobs = cursor.fetchall()
    if jobs:
        accessions = [str(k[1]) for k in jobs]
        files = [str(k[2]) for k in jobs]
        return {taxon: accessions}, {taxon: files}
    else:
        return None, None


def retrieve_taxa_analyzed(con, taxon):
    query = "SELECT * from genomes where taxon =?"
    cursor = con.cursor()
    cursor.execute(query, (taxon,))
    jobs = cursor.fetchall()
    if jobs:
        accessions = [str(k[1]) for k in jobs]
        return {taxon: accessions}
    else:
        return None


def retrieve_taxa_analyzed(con, taxon):
    query = "SELECT * from genomes where taxon =?"
    cursor = con.cursor()
    cursor.execute(query, (taxon,))
    jobs = cursor.fetchall()
    if jobs:
        accessions = [str(k[1]) for k in jobs]
        return {taxon: accessions}
    else:
        return None


def retrieve_results_done(con, taxon):
    query = "SELECT * from results where taxon =?"
    cursor = con.cursor()
    cursor.execute(query, (taxon,))
    jobs = cursor.fetchall()
    if jobs:
        return True
    else:
        return False


def check_if_updates(genomes_done, acc2assm):
    if set(next(iter(acc2assm.keys()))) == set(next(iter(genomes_done.values()))):
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
