"""
DuckDB connection management and utility functions.
"""

import os
import logging
from typing import Dict, List, Optional, Any, Union, Tuple
from pathlib import Path

import duckdb
import numpy as np

logger = logging.getLogger(__name__)


class DuckDBManager:
    """
    Manager for DuckDB connections and operations.
    """

    def __init__(self, db_path: str, read_only: bool = False):
        """
        Initialize DuckDB connection.

        Args:
            db_path: Path to DuckDB database file
            read_only: Whether to open in read-only mode
        """
        self.db_path = db_path
        self.read_only = read_only
        self.conn = None

        # Connect to database
        self._connect()

    def _connect(self):
        """Connect to DuckDB database."""
        try:
            self.conn = duckdb.connect(self.db_path, read_only=self.read_only)
            logger.debug(f"Connected to DuckDB at {self.db_path}")
        except Exception as e:
            logger.error(f"Failed to connect to DuckDB at {self.db_path}: {e}")
            raise

    def close(self):
        """Close DuckDB connection."""
        if self.conn:
            self.conn.close()
            self.conn = None
            logger.debug(f"Closed connection to DuckDB at {self.db_path}")

    def __enter__(self):
        """Context manager entry."""
        if not self.conn:
            self._connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()

    def execute(self, query: str, parameters: Optional[Dict[str, Any]] = None) -> Any:
        """
        Execute a SQL query.

        Args:
            query: SQL query string
            parameters: Query parameters

        Returns:
            Query result
        """
        try:
            if parameters:
                result = self.conn.execute(query, parameters)
            else:
                result = self.conn.execute(query)
            return result
        except Exception as e:
            logger.error(f"Failed to execute query: {e}")
            logger.debug(f"Query: {query}")
            logger.debug(f"Parameters: {parameters}")
            raise

    def fetch_all(
        self, query: str, parameters: Optional[Dict[str, Any]] = None
    ) -> List[Tuple]:
        """
        Execute a query and fetch all results.

        Args:
            query: SQL query string
            parameters: Query parameters

        Returns:
            List of result tuples
        """
        result = self.execute(query, parameters)
        return result.fetchall()

    def fetch_one(
        self, query: str, parameters: Optional[Dict[str, Any]] = None
    ) -> Tuple:
        """
        Execute a query and fetch one result.

        Args:
            query: SQL query string
            parameters: Query parameters

        Returns:
            Result tuple
        """
        result = self.execute(query, parameters)
        return result.fetchone()

    def create_table(
        self,
        table_name: str,
        columns: Dict[str, str],
        primary_key: Optional[str] = None,
    ):
        """
        Create a table if it doesn't exist.

        Args:
            table_name: Name of the table
            columns: Dictionary of column names and types
            primary_key: Primary key column name
        """
        # Construct column definitions
        col_defs = []
        for col_name, col_type in columns.items():
            if primary_key and col_name == primary_key:
                col_defs.append(f"{col_name} {col_type} PRIMARY KEY")
            else:
                col_defs.append(f"{col_name} {col_type}")

        # Create table
        query = f"CREATE TABLE IF NOT EXISTS {table_name} ({', '.join(col_defs)})"
        self.execute(query)
        logger.debug(f"Created table {table_name} if it didn't exist")

    def create_index(
        self, table_name: str, column_name: str, index_name: Optional[str] = None
    ):
        """
        Create an index on a table column.

        Args:
            table_name: Name of the table
            column_name: Name of the column to index
            index_name: Name of the index (optional)
        """
        if not index_name:
            index_name = f"idx_{table_name}_{column_name}"

        query = (
            f"CREATE INDEX IF NOT EXISTS {index_name} ON {table_name} ({column_name})"
        )
        self.execute(query)
        logger.debug(f"Created index {index_name} on {table_name}.{column_name}")

    def insert(self, table_name: str, data: Dict[str, List[Any]]):
        """
        Insert data into a table.

        Args:
            table_name: Name of the table
            data: Dictionary of column names and values
        """
        # Check if data is empty
        if not data or not data[list(data.keys())[0]]:
            logger.debug(f"No data to insert into {table_name}")
            return

        # Get column names
        columns = list(data.keys())

        # Get number of rows
        n_rows = len(data[columns[0]])

        # Prepare parameter names for query
        param_names = [f"${i+1}" for i in range(len(columns))]

        # Construct query
        query = f"INSERT INTO {table_name} ({', '.join(columns)}) VALUES ({', '.join(param_names)})"

        # Insert rows one by one
        # TODO: Use bulk insert for better performance
        inserted = 0

        for i in range(n_rows):
            params = [data[col][i] for col in columns]
            try:
                self.execute(query, params)
                inserted += 1
            except Exception as e:
                logger.warning(f"Failed to insert row {i} into {table_name}: {e}")

        logger.debug(f"Inserted {inserted} rows into {table_name}")

    def bulk_insert(
        self, table_name: str, data: Dict[str, List[Any]], batch_size: int = 1000
    ):
        """
        Insert data into a table in batches.

        Args:
            table_name: Name of the table
            data: Dictionary of column names and values
            batch_size: Number of rows per batch
        """
        # Check if data is empty
        if not data or not data[list(data.keys())[0]]:
            logger.debug(f"No data to insert into {table_name}")
            return

        # Get column names
        columns = list(data.keys())

        # Get number of rows
        n_rows = len(data[columns[0]])

        # Calculate number of batches
        n_batches = (n_rows + batch_size - 1) // batch_size

        inserted = 0

        for batch_idx in range(n_batches):
            # Get start and end indices for this batch
            start_idx = batch_idx * batch_size
            end_idx = min((batch_idx + 1) * batch_size, n_rows)

            # Prepare parameter names for query
            param_names = [f"${i+1}" for i in range(len(columns))]

            # Construct query for this batch
            query = f"INSERT INTO {table_name} ({', '.join(columns)}) VALUES "

            # Add value placeholders for each row in the batch
            values_clauses = []
            for _ in range(start_idx, end_idx):
                values_clauses.append(f"({', '.join(param_names)})")

            query += ", ".join(values_clauses)

            # Flatten parameters for all rows in the batch
            params = []
            for i in range(start_idx, end_idx):
                for col in columns:
                    params.append(data[col][i])

            try:
                self.execute(query, params)
                inserted += end_idx - start_idx
            except Exception as e:
                logger.warning(
                    f"Failed to insert batch {batch_idx} into {table_name}: {e}"
                )

        logger.debug(
            f"Inserted {inserted} rows into {table_name} in {n_batches} batches"
        )

    def table_exists(self, table_name: str) -> bool:
        """
        Check if a table exists.

        Args:
            table_name: Name of the table

        Returns:
            True if the table exists, False otherwise
        """
        query = f"SELECT count(*) FROM information_schema.tables WHERE table_name = '{table_name}'"
        result = self.fetch_one(query)
        return result[0] > 0

    def get_tables(self) -> List[str]:
        """
        Get a list of all tables in the database.

        Returns:
            List of table names
        """
        query = "SELECT table_name FROM information_schema.tables"
        result = self.fetch_all(query)
        return [row[0] for row in result]

    def drop_table(self, table_name: str):
        """
        Drop a table.

        Args:
            table_name: Name of the table
        """
        query = f"DROP TABLE IF EXISTS {table_name}"
        self.execute(query)
        logger.debug(f"Dropped table {table_name}")

    def query_to_dict(
        self, query: str, parameters: Optional[Dict[str, Any]] = None
    ) -> Dict[str, List[Any]]:
        """
        Execute a query and return results as a dictionary.

        Args:
            query: SQL query string
            parameters: Query parameters

        Returns:
            Dictionary of column names and values
        """
        result = self.execute(query, parameters)

        # Get column names
        column_names = [desc[0] for desc in result.description]

        # Initialize result dictionary
        result_dict = {col: [] for col in column_names}

        # Fill result dictionary
        for row in result.fetchall():
            for i, col in enumerate(column_names):
                result_dict[col].append(row[i])

        return result_dict

    def begin_transaction(self):
        """Begin a transaction."""
        self.execute("BEGIN TRANSACTION")

    def commit(self):
        """Commit the current transaction."""
        self.execute("COMMIT")

    def rollback(self):
        """Rollback the current transaction."""
        self.execute("ROLLBACK")

    def vacuum(self):
        """Vacuum the database to optimize storage."""
        self.execute("VACUUM")
        logger.debug("Vacuumed database")
