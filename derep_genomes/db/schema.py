"""
Database schema definitions for the genome dereplication database.
"""

import logging
from typing import Dict, List, Optional, Any, Union

from derep_genomes.db.duckdb_manager import DuckDBManager

logger = logging.getLogger(__name__)


# Schema definitions for all tables
SCHEMA = {
    "taxa": {
        "columns": {"taxonomy": "VARCHAR", "count": "INTEGER"},
        "primary_key": "taxonomy",
        "indices": ["taxonomy"],
    },
    "genomes": {
        "columns": {
            "taxonomy": "VARCHAR",
            "accession": "VARCHAR",
            "file": "VARCHAR",
            "genome_length": "INTEGER",
        },
        "primary_key": "accession",
        "indices": ["taxonomy", "accession"],
    },
    "genomes_derep": {
        "columns": {
            "accession": "VARCHAR",
            "representative": "BOOLEAN",
            "derep": "BOOLEAN",
        },
        "primary_key": "accession",
        "indices": ["accession", "representative", "derep"],
    },
    "results": {
        "columns": {
            "taxonomy": "VARCHAR",
            "weight": "FLOAT",
            "communities": "INTEGER",
            "n_genomes": "INTEGER",
            "n_genomes_derep": "INTEGER",
            "method": "VARCHAR",
        },
        "primary_key": "taxonomy",
        "indices": ["taxonomy", "method"],
    },
    "stats": {
        "columns": {
            "taxonomy": "VARCHAR",
            "representative": "VARCHAR",
            "n_nodes": "INTEGER",
            "n_nodes_selected": "INTEGER",
            "n_nodes_discarded": "INTEGER",
            "graph_avg_weight": "FLOAT",
            "graph_sd_weight": "FLOAT",
            "graph_avg_weight_raw": "FLOAT",
            "graph_sd_weight_raw": "FLOAT",
            "subgraph_selected_avg_weight": "FLOAT",
            "subgraph_selected_sd_weight": "FLOAT",
            "subgraph_selected_avg_weight_raw": "FLOAT",
            "subgraph_selected_sd_weight_raw": "FLOAT",
            "subgraph_discarded_avg_weight": "FLOAT",
            "subgraph_discarded_sd_weight": "FLOAT",
            "subgraph_discarded_avg_weight_raw": "FLOAT",
            "subgraph_discarded_sd_weight_raw": "FLOAT",
        },
        "primary_key": "representative",
        "indices": ["taxonomy", "representative"],
    },
    "jobs_done": {
        "columns": {"taxonomy": "VARCHAR", "accession": "VARCHAR", "file": "VARCHAR"},
        "primary_key": "accession",
        "indices": ["taxonomy", "accession"],
    },
    "jobs_failed": {
        "columns": {
            "taxonomy": "VARCHAR",
            "accession": "VARCHAR",
            "file": "VARCHAR",
            "reason": "VARCHAR",
        },
        "primary_key": "accession",
        "indices": ["taxonomy", "accession"],
    },
}


def create_schema(db_manager: DuckDBManager) -> None:
    """
    Create database schema (tables and indices).

    Args:
        db_manager: DuckDB connection manager
    """
    # Create tables from schema
    for table_name, table_schema in SCHEMA.items():
        # Create table
        db_manager.create_table(
            table_name=table_name,
            columns=table_schema["columns"],
            primary_key=table_schema.get("primary_key"),
        )

        # Create indices
        for column in table_schema.get("indices", []):
            db_manager.create_index(
                table_name=table_name,
                column_name=column,
                index_name=f"idx_{table_name}_{column}",
            )

    logger.info("Database schema created successfully")


def check_schema(db_manager: DuckDBManager) -> bool:
    """
    Check if the database schema exists and is valid.

    Args:
        db_manager: DuckDB connection manager

    Returns:
        True if schema exists and is valid, False otherwise
    """
    # Get list of existing tables
    existing_tables = db_manager.get_tables()

    # Check if all tables exist
    for table_name in SCHEMA.keys():
        if table_name.lower() not in [t.lower() for t in existing_tables]:
            logger.warning(f"Table {table_name} does not exist")
            return False

    # Schema is valid
    return True


def reset_schema(db_manager: DuckDBManager) -> None:
    """
    Reset database schema by dropping all tables and recreating them.

    Args:
        db_manager: DuckDB connection manager
    """
    # Drop all tables
    for table_name in SCHEMA.keys():
        db_manager.drop_table(table_name)

    # Recreate schema
    create_schema(db_manager)

    logger.info("Database schema reset successfully")
