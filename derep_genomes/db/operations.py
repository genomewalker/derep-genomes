"""
Database operations for genome dereplication.
"""

import logging
from typing import Dict, List, Optional, Any, Union, Set, Tuple

from derep_genomes.db.duckdb_manager import DuckDBManager

logger = logging.getLogger(__name__)


class DBOperations:
    """Database operations for genome dereplication."""

    def __init__(self, db_manager: DuckDBManager):
        """
        Initialize database operations.

        Args:
            db_manager: DuckDB connection manager
        """
        self.db = db_manager

    def insert_taxon(self, taxonomy: str, count: int) -> None:
        """
        Insert a taxon into the database.

        Args:
            taxonomy: Taxon name
            count: Number of genomes in the taxon
        """
        # Check if taxon already exists
        query = "SELECT COUNT(*) FROM taxa WHERE taxonomy = ?"
        result = self.db.fetch_one(query, [taxonomy])

        if result[0] > 0:
            # Taxon already exists, update count
            query = "UPDATE taxa SET count = ? WHERE taxonomy = ?"
            self.db.execute(query, [count, taxonomy])
            logger.debug(f"Updated taxon {taxonomy} with count {count}")
        else:
            # Insert new taxon
            query = "INSERT INTO taxa (taxonomy, count) VALUES (?, ?)"
            self.db.execute(query, [taxonomy, count])
            logger.debug(f"Inserted taxon {taxonomy} with count {count}")

    def insert_genomes(self, genomes: List[Dict[str, Any]]) -> None:
        """
        Insert genomes into the database.

        Args:
            genomes: List of genome dictionaries with keys:
                - taxonomy
                - accession
                - file
                - genome_length
        """
        if not genomes:
            logger.debug("No genomes to insert")
            return

        # Convert list of dictionaries to dictionary of lists
        data = {"taxonomy": [], "accession": [], "file": [], "genome_length": []}

        for genome in genomes:
            data["taxonomy"].append(genome["taxonomy"])
            data["accession"].append(genome["accession"])
            data["file"].append(genome["file"])
            data["genome_length"].append(genome.get("genome_length", 0))

        # Bulk insert
        self.db.bulk_insert("genomes", data)

    def insert_derep_genomes(self, derep_genomes: List[Dict[str, Any]]) -> None:
        """
        Insert dereplicated genomes into the database.

        Args:
            derep_genomes: List of dereplicated genome dictionaries with keys:
                - accession
                - representative
                - derep
        """
        if not derep_genomes:
            logger.debug("No dereplicated genomes to insert")
            return

        # Convert list of dictionaries to dictionary of lists
        data = {"accession": [], "representative": [], "derep": []}

        for genome in derep_genomes:
            data["accession"].append(genome["accession"])
            # Convert to boolean
            data["representative"].append(bool(genome["representative"]))
            data["derep"].append(bool(genome["derep"]))

        # Bulk insert
        self.db.bulk_insert("genomes_derep", data)

    def insert_results(self, results: Dict[str, Any]) -> None:
        """
        Insert dereplication results into the database.

        Args:
            results: Results dictionary with keys:
                - taxonomy
                - weight
                - communities
                - n_genomes
                - n_genomes_derep
                - method
        """
        # Check if results already exist for this taxon
        query = "SELECT COUNT(*) FROM results WHERE taxonomy = ? AND method = ?"
        result = self.db.fetch_one(query, [results["taxonomy"], results["method"]])

        if result[0] > 0:
            # Results already exist, update
            query = """
                UPDATE results 
                SET weight = ?, communities = ?, n_genomes = ?, n_genomes_derep = ?
                WHERE taxonomy = ? AND method = ?
            """
            self.db.execute(
                query,
                [
                    results["weight"],
                    results["communities"],
                    results["n_genomes"],
                    results["n_genomes_derep"],
                    results["taxonomy"],
                    results["method"],
                ],
            )
            logger.debug(
                f"Updated results for taxon {results['taxonomy']} using method {results['method']}"
            )
        else:
            # Insert new results
            query = """
                INSERT INTO results 
                (taxonomy, weight, communities, n_genomes, n_genomes_derep, method)
                VALUES (?, ?, ?, ?, ?, ?)
            """
            self.db.execute(
                query,
                [
                    results["taxonomy"],
                    results["weight"],
                    results["communities"],
                    results["n_genomes"],
                    results["n_genomes_derep"],
                    results["method"],
                ],
            )
            logger.debug(
                f"Inserted results for taxon {results['taxonomy']} using method {results['method']}"
            )

    def insert_stats(self, stats: List[Dict[str, Any]]) -> None:
        """
        Insert statistics into the database.

        Args:
            stats: List of statistics dictionaries
        """
        if not stats:
            logger.debug("No statistics to insert")
            return

        # Convert list of dictionaries to dictionary of lists
        data = {
            "taxonomy": [],
            "representative": [],
            "n_nodes": [],
            "n_nodes_selected": [],
            "n_nodes_discarded": [],
            "graph_avg_weight": [],
            "graph_sd_weight": [],
            "graph_avg_weight_raw": [],
            "graph_sd_weight_raw": [],
            "subgraph_selected_avg_weight": [],
            "subgraph_selected_sd_weight": [],
            "subgraph_selected_avg_weight_raw": [],
            "subgraph_selected_sd_weight_raw": [],
            "subgraph_discarded_avg_weight": [],
            "subgraph_discarded_sd_weight": [],
            "subgraph_discarded_avg_weight_raw": [],
            "subgraph_discarded_sd_weight_raw": [],
        }

        for stat in stats:
            for key in data:
                data[key].append(stat.get(key, None))

        # Bulk insert
        self.db.bulk_insert("stats", data)

    def insert_job_done(self, job: Dict[str, Any]) -> None:
        """
        Insert a completed job into the database.

        Args:
            job: Job dictionary with keys:
                - taxonomy
                - accession
                - file
        """
        query = "INSERT INTO jobs_done (taxonomy, accession, file) VALUES (?, ?, ?)"
        self.db.execute(query, [job["taxonomy"], job["accession"], job["file"]])
        logger.debug(f"Inserted completed job for {job['accession']}")

    def insert_job_failed(self, job: Dict[str, Any]) -> None:
        """
        Insert a failed job into the database.

        Args:
            job: Job dictionary with keys:
                - taxonomy
                - accession
                - file
                - reason
        """
        query = "INSERT INTO jobs_failed (taxonomy, accession, file, reason) VALUES (?, ?, ?, ?)"
        self.db.execute(
            query, [job["taxonomy"], job["accession"], job["file"], job["reason"]]
        )
        logger.debug(
            f"Inserted failed job for {job['accession']} with reason: {job['reason']}"
        )

    def get_taxon_genomes(self, taxonomy: str) -> List[Dict[str, Any]]:
        """
        Get all genomes for a taxon.

        Args:
            taxonomy: Taxon name

        Returns:
            List of genome dictionaries
        """
        query = "SELECT * FROM genomes WHERE taxonomy = ?"
        result = self.db.fetch_all(query, [taxonomy])

        genomes = []
        for row in result:
            genomes.append(
                {
                    "taxonomy": row[0],
                    "accession": row[1],
                    "file": row[2],
                    "genome_length": row[3],
                }
            )

        return genomes

    def get_all_taxa(self) -> List[Dict[str, Any]]:
        """
        Get all taxa in the database.

        Returns:
            List of taxon dictionaries
        """
        query = "SELECT * FROM taxa"
        result = self.db.fetch_all(query)

        taxa = []
        for row in result:
            taxa.append({"taxonomy": row[0], "count": row[1]})

        return taxa

    def get_derep_genomes(
        self, representatives_only: bool = False
    ) -> List[Dict[str, Any]]:
        """
        Get all dereplicated genomes.

        Args:
            representatives_only: Whether to get only representative genomes

        Returns:
            List of dereplicated genome dictionaries
        """
        if representatives_only:
            query = "SELECT * FROM genomes_derep WHERE representative = TRUE"
        else:
            query = "SELECT * FROM genomes_derep"

        result = self.db.fetch_all(query)

        genomes = []
        for row in result:
            genomes.append(
                {"accession": row[0], "representative": row[1], "derep": row[2]}
            )

        return genomes

    def get_all_jobs_done(self) -> List[Dict[str, Any]]:
        """
        Get all completed jobs.

        Returns:
            List of completed job dictionaries
        """
        query = "SELECT * FROM jobs_done"
        result = self.db.fetch_all(query)

        jobs = []
        for row in result:
            jobs.append({"taxonomy": row[0], "accession": row[1], "file": row[2]})

        return jobs

    def get_all_jobs_failed(self) -> List[Dict[str, Any]]:
        """
        Get all failed jobs.

        Returns:
            List of failed job dictionaries
        """
        query = "SELECT * FROM jobs_failed"
        result = self.db.fetch_all(query)

        jobs = []
        for row in result:
            jobs.append(
                {
                    "taxonomy": row[0],
                    "accession": row[1],
                    "file": row[2],
                    "reason": row[3],
                }
            )

        return jobs

    def has_taxon_been_processed(self, taxonomy: str) -> bool:
        """
        Check if a taxon has been processed.

        Args:
            taxonomy: Taxon name

        Returns:
            True if the taxon has been processed, False otherwise
        """
        query = "SELECT COUNT(*) FROM results WHERE taxonomy = ?"
        result = self.db.fetch_one(query, [taxonomy])

        return result[0] > 0

    def delete_taxon(self, taxonomy: str) -> None:
        """
        Delete a taxon and all associated data from the database.

        Args:
            taxonomy: Taxon name
        """
        # Start transaction
        self.db.begin_transaction()

        try:
            # Get all accessions for this taxon
            query = "SELECT accession FROM genomes WHERE taxonomy = ?"
            result = self.db.fetch_all(query, [taxonomy])
            accessions = [row[0] for row in result]

            # Create placeholder string for IN clause
            placeholders = ", ".join(["?"] * len(accessions)) if accessions else "''"

            # Delete from jobs_done
            query = f"DELETE FROM jobs_done WHERE taxonomy = ?"
            self.db.execute(query, [taxonomy])

            # Delete from jobs_failed
            query = f"DELETE FROM jobs_failed WHERE taxonomy = ?"
            self.db.execute(query, [taxonomy])

            # Delete from stats
            query = f"DELETE FROM stats WHERE taxonomy = ?"
            self.db.execute(query, [taxonomy])

            # Delete from results
            query = f"DELETE FROM results WHERE taxonomy = ?"
            self.db.execute(query, [taxonomy])

            # Delete from genomes_derep
            if accessions:
                query = f"DELETE FROM genomes_derep WHERE accession IN ({placeholders})"
                self.db.execute(query, accessions)

            # Delete from genomes
            query = f"DELETE FROM genomes WHERE taxonomy = ?"
            self.db.execute(query, [taxonomy])

            # Delete from taxa
            query = f"DELETE FROM taxa WHERE taxonomy = ?"
            self.db.execute(query, [taxonomy])

            # Commit transaction
            self.db.commit()

            logger.info(f"Deleted taxon {taxonomy} and all associated data")

        except Exception as e:
            # Rollback transaction on error
            self.db.rollback()
            logger.error(f"Failed to delete taxon {taxonomy}: {e}")
            raise
