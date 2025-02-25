"""
Main entry point for the genome dereplication tool.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

import os
import sys
import shutil
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Set, Tuple

from derep_genomes.config import Config, parse_args
from derep_genomes.db.duckdb_manager import DuckDBManager
from derep_genomes.db.schema import create_schema, check_schema
from derep_genomes.db.operations import DBOperations
from derep_genomes.utils.logging import (
    setup_logging,
    Timer,
    log_memory_usage,
    is_debug_enabled,
)
from derep_genomes.utils.parallel import ParallelTaskManager
from derep_genomes.utils.helpers import split_into_chunks
from derep_genomes.core.similarity import SimilarityMeasure

# Set up logger
logger = logging.getLogger("derep_genomes")


def setup_run_directories(config: Config) -> Dict[str, Path]:
    """
    Create and setup directories for the current run.

    Args:
        config: Configuration object

    Returns:
        Dictionary of directory paths
    """
    # Create base run directory with timestamp
    base_dir = (
        Path(config.out_dir if config.out_dir else ".")
        / f"derep_run_{config.timestamp}"
    )

    # Create directory structure
    dirs = {
        "base": base_dir,
        "temp": base_dir / "temp",
        "results": base_dir / "results",
        "logs": base_dir / "logs",
        "reps": base_dir / "representatives" if config.copy else None,
    }

    # Create directories
    for name, path in dirs.items():
        if path:
            os.makedirs(path, exist_ok=True)
            logger.debug(f"Created directory: {path}")

    # Set output paths in config
    config.results_dir = str(dirs["results"])
    config.temp_dir = str(dirs["temp"])
    if config.copy:
        config.rep_out_dir = str(dirs["reps"])

    # Update log file location if not explicitly set
    if not config.log_file:
        config.log_file = str(dirs["logs"] / f"derep_{config.timestamp}.log")
        # Re-setup logging with file
        from derep_genomes.utils.logging import setup_logging

        log_level = logging.DEBUG if config.debug else logging.INFO
        setup_logging(log_level=log_level, log_file=config.log_file)

    return dirs


def sanitize_for_filename(text: str) -> str:
    """
    Sanitize a string to make it safe for use in a filename.

    Args:
        text: Input string

    Returns:
        Sanitized string
    """
    # Replace problematic characters
    replacements = {
        ";": "_",
        " ": "_",
        ":": "_",
        "/": "_",
        "\\": "_",
        "*": "_",
        "?": "_",
        '"': "_",
        "<": "_",
        ">": "_",
        "|": "_",
    }

    for char, replacement in replacements.items():
        text = text.replace(char, replacement)

    return text


def load_taxa_from_file(file_path: str) -> List[str]:
    """
    Load taxa from a file.

    Args:
        file_path: Path to taxa file

    Returns:
        List of taxa
    """
    taxa = []

    try:
        with open(file_path, "r") as f:
            for line in f:
                taxa.append(line.strip())
    except Exception as e:
        logger.error(f"Error loading taxa from {file_path}: {e}")
        sys.exit(1)

    return taxa


def load_genome_data(file_path: str) -> List[Dict[str, Any]]:
    """
    Load genome data from a TSV file.

    Args:
        file_path: Path to genome data file

    Returns:
        List of genome dictionaries
    """
    genomes = []

    try:
        with open(file_path, "r") as f:
            # Skip header if present
            header = f.readline().strip().split("\t")

            # Check if header contains expected columns
            if len(header) < 2 or "accession" not in header and "file" not in header:
                # No header, reset file pointer
                f.seek(0)
                header = ["accession", "taxonomy", "file"]

            # Get column indices
            acc_idx = header.index("accession") if "accession" in header else 0
            tax_idx = header.index("taxonomy") if "taxonomy" in header else 1
            file_idx = header.index("file") if "file" in header else 2

            # Read data
            for line in f:
                fields = line.strip().split("\t")

                if len(fields) <= max(acc_idx, tax_idx, file_idx):
                    logger.warning(f"Skipping invalid line: {line.strip()}")
                    continue

                genome = {
                    "accession": fields[acc_idx],
                    "taxonomy": fields[tax_idx],
                    "file": fields[file_idx],
                }

                genomes.append(genome)
    except Exception as e:
        logger.error(f"Error loading genome data from {file_path}: {e}")
        sys.exit(1)

    return genomes


def filter_by_taxa(
    genomes: List[Dict[str, Any]], taxa: List[str]
) -> List[Dict[str, Any]]:
    """
    Filter genomes by taxa.

    Args:
        genomes: List of genome dictionaries
        taxa: List of taxa to keep

    Returns:
        Filtered list of genome dictionaries
    """
    return [genome for genome in genomes if genome["taxonomy"] in taxa]


def check_files_exist(genomes: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Check if genome files exist and filter out those that don't.

    Args:
        genomes: List of genome dictionaries

    Returns:
        Filtered list of genome dictionaries
    """
    valid_genomes = []

    for genome in genomes:
        if os.path.isfile(genome["file"]):
            valid_genomes.append(genome)
        else:
            logger.warning(
                f"File not found: {genome['file']} for accession {genome['accession']}"
            )

    return valid_genomes


def group_by_taxonomy(genomes: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    """
    Group genomes by taxonomy.

    Args:
        genomes: List of genome dictionaries

    Returns:
        Dictionary of taxonomy to list of genome dictionaries
    """
    taxonomy_groups = {}

    for genome in genomes:
        taxonomy = genome["taxonomy"]
        if taxonomy not in taxonomy_groups:
            taxonomy_groups[taxonomy] = []
        taxonomy_groups[taxonomy].append(genome)

    return taxonomy_groups


def get_similarity_measure(config: Config, temp_dir: Path) -> SimilarityMeasure:
    """
    Get appropriate similarity measure based on configuration.

    Args:
        config: Configuration object
        temp_dir: Temporary directory

    Returns:
        Similarity measure instance
    """
    if config.similarity_method == "skani":
        from derep_genomes.core.skani import SKANI

        return SKANI(
            threads=config.threads,
            threshold=config.threshold,
            temp_dir=temp_dir,
            skani_k=config.skani_k,
            skani_s=config.skani_s,
            skani_min_conf=config.skani_min_conf,
            skani_target_ani=config.skani_target_ani,
            skani_min_kmers=config.skani_min_kmers,
        )
    else:  # default to FastANI
        from derep_genomes.core.ani import FastANI

        return FastANI(
            threads=config.threads,
            threshold=config.threshold,
            temp_dir=temp_dir,
            ani_fraglen=config.ani_fraglen,
            ani_fraglen_fraction=config.ani_fraglen_fraction,
            ani_min_frac=config.ani_min_frac,
            ani_threshold=config.ani_threshold,
            min_genome_size=config.min_genome_size,
        )


def get_precluster_similarity_measure(
    config: Config, temp_dir: Path
) -> SimilarityMeasure:
    """
    Get appropriate pre-clustering similarity measure based on configuration.

    Args:
        config: Configuration object
        temp_dir: Temporary directory

    Returns:
        Similarity measure instance
    """
    if config.dashing:
        from derep_genomes.core.dashing import Dashing

        return Dashing(
            threads=config.threads,
            threshold=config.threshold,
            temp_dir=temp_dir,
            dashing_k=config.dashing_k,
            dashing_s=config.dashing_s,
            dashing_threshold=config.dashing_threshold,
        )
    else:  # default to Mash
        from derep_genomes.core.mash import Mash

        return Mash(
            threads=config.threads,
            threshold=config.threshold,
            temp_dir=temp_dir,
            mash_k=config.mash_k,
            mash_s=config.mash_s,
            mash_threshold=config.mash_threshold,
        )


def process_taxon(
    taxon: str, genomes: List[Dict[str, Any]], config: Config, db_ops: DBOperations
) -> Dict[str, Any]:
    """
    Process a single taxon.

    Args:
        taxon: Taxonomy name
        genomes: List of genome dictionaries for this taxon
        config: Configuration object
        db_ops: Database operations

    Returns:
        Processing results
    """
    logger.info(f"Processing taxon: {taxon} with {len(genomes)} genomes")

    # Skip if already processed
    if db_ops.has_taxon_been_processed(taxon):
        logger.info(f"Taxon {taxon} already processed, skipping")
        return {"taxonomy": taxon, "status": "skipped", "n_genomes": len(genomes)}

    # If singleton taxon, mark all genomes as representatives
    if len(genomes) == 1:
        logger.info(f"Taxon {taxon} has only one genome, marking as representative")

        # Insert into database
        db_ops.insert_taxon(taxon, 1)
        db_ops.insert_genomes(genomes)

        derep_genome = {
            "accession": genomes[0]["accession"],
            "representative": True,
            "derep": True,
        }

        db_ops.insert_derep_genomes([derep_genome])

        result = {
            "taxonomy": taxon,
            "weight": None,
            "communities": 1,
            "n_genomes": 1,
            "n_genomes_derep": 1,
            "method": config.similarity_method,
        }

        db_ops.insert_results(result)

        for genome in genomes:
            db_ops.insert_job_done(
                {
                    "taxonomy": taxon,
                    "accession": genome["accession"],
                    "file": genome["file"],
                }
            )

        return {
            "taxonomy": taxon,
            "status": "singleton",
            "n_genomes": 1,
            "representatives": [genomes[0]["file"]],
        }

    # Create temporary directory with sanitized taxon name
    safe_taxon = sanitize_for_filename(taxon)
    temp_dir = Path(config.tmp_dir) / f"derep-{safe_taxon}-{config.timestamp}"
    os.makedirs(temp_dir, exist_ok=True)

    try:
        # Extract assembly file paths
        assembly_files = [genome["file"] for genome in genomes]

        # Prefetch genome lengths for all assembly files at the beginning
        from derep_genomes.utils.helpers import genome_length_cache

        genome_length_cache.prefetch(assembly_files)

        # First stage: Pre-clustering with Mash or Dashing
        logger.info(
            f"Performing pre-clustering with {'Dashing' if config.dashing else 'Mash'}"
        )

        # Get pre-clustering similarity measure
        precluster_measure = get_precluster_similarity_measure(config, temp_dir)

        with Timer(f"Pre-clustering for {taxon}"):
            precluster_data = precluster_measure.compute_pairwise(assembly_files)
            filtered_precluster_data = precluster_measure.filter_by_similarity(
                precluster_data
            )

            # Check for missing assemblies (too divergent)
            missing_assemblies = precluster_measure.check_missing_assemblies(
                filtered_precluster_data, assembly_files
            )

            # Create graph for pre-clustered data
            precluster_graph = precluster_measure.create_similarity_graph(
                filtered_precluster_data
            )

        # Find connected components in pre-clustering graph
        from networkx import connected_components
        import networkx as nx

        # Group assemblies by connected component
        clusters = []
        for component in connected_components(precluster_graph):
            if len(component) > 0:
                clusters.append(list(component))

        logger.info(f"Pre-clustering identified {len(clusters)} clusters")

        # Second stage: Detailed similarity analysis within each cluster
        final_representatives = []

        for i, cluster in enumerate(clusters):
            if len(cluster) == 1:
                # Single assembly in cluster, add as representative
                final_representatives.extend(cluster)
                continue

            logger.info(
                f"Processing cluster {i+1}/{len(clusters)} with {len(cluster)} genomes"
            )

            # Get appropriate similarity measure for detailed analysis
            similarity_measure = get_similarity_measure(config, temp_dir)

            with Timer(f"Similarity calculation for cluster {i+1}"):
                similarity_data = similarity_measure.compute_pairwise(cluster)
                filtered_data = similarity_measure.filter_by_similarity(similarity_data)

                # Create graph and find representative
                from derep_genomes.core.graph import find_clusters

                graph = similarity_measure.create_similarity_graph(filtered_data)
                cluster_results = find_clusters(graph, threshold=config.threshold)

                # Add representatives from this cluster
                final_representatives.extend(cluster_results["representatives"])

        # Add missing assemblies (too divergent) as representatives
        final_representatives.extend(missing_assemblies)

        # Remove duplicates
        final_representatives = list(set(final_representatives))

        # Create results
        logger.info(
            f"Found {len(final_representatives)} representatives out of {len(genomes)} genomes"
        )

        # Insert into database
        db_ops.insert_taxon(taxon, len(genomes))
        db_ops.insert_genomes(genomes)

        # Insert dereplicated genomes
        derep_genomes = []
        for genome in genomes:
            is_rep = genome["file"] in final_representatives

            derep_genomes.append(
                {
                    "accession": genome["accession"],
                    "representative": is_rep,
                    "derep": True,  # All are part of the dereplicated set
                }
            )

        db_ops.insert_derep_genomes(derep_genomes)

        # Insert results
        result = {
            "taxonomy": taxon,
            "weight": None,  # Not applicable with two-stage approach
            "communities": len(clusters),
            "n_genomes": len(genomes),
            "n_genomes_derep": len(final_representatives),
            "method": f"{'dashing' if config.dashing else 'mash'}+{config.similarity_method}",
        }

        db_ops.insert_results(result)

        # Insert jobs done
        for genome in genomes:
            db_ops.insert_job_done(
                {
                    "taxonomy": taxon,
                    "accession": genome["accession"],
                    "file": genome["file"],
                }
            )

        # Get memory usage statistics
        if is_debug_enabled():
            log_memory_usage(f"Memory usage after processing {taxon}")
            logger.debug(
                f"Genome length cache stats: {genome_length_cache.get_stats()}"
            )

        # Clean up
        logger.debug(f"Cleaning up temporary directory: {temp_dir}")
        shutil.rmtree(temp_dir, ignore_errors=True)

        return {
            "taxonomy": taxon,
            "status": "success",
            "n_genomes": len(genomes),
            "n_representatives": len(final_representatives),
            "representatives": final_representatives,
        }

    except Exception as e:
        logger.error(f"Error processing taxon {taxon}: {e}")

        # Mark job as failed
        for genome in genomes:
            db_ops.insert_job_failed(
                {
                    "taxonomy": taxon,
                    "accession": genome["accession"],
                    "file": genome["file"],
                    "reason": str(e),
                }
            )

        # Clean up
        logger.debug(f"Cleaning up temporary directory: {temp_dir}")
        shutil.rmtree(temp_dir, ignore_errors=True)

        return {
            "taxonomy": taxon,
            "status": "failed",
            "n_genomes": len(genomes),
            "error": str(e),
        }


def copy_representative_files(
    results: List[Dict[str, Any]], genomes: List[Dict[str, Any]], out_dir: str
) -> None:
    """
    Copy representative genome files to output directory.

    Args:
        results: List of processing results
        genomes: List of all genome dictionaries
        out_dir: Output directory
    """
    # Create output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)

    # Create lookup for genomes by file path
    genome_lookup = {genome["file"]: genome for genome in genomes}

    # Collect all representative files
    rep_files = []
    for result in results:
        if result.get("status") == "success" or result.get("status") == "singleton":
            rep_files.extend(result.get("representatives", []))

    # Copy files
    logger.info(f"Copying {len(rep_files)} representative genomes to {out_dir}")

    with Timer("Copying representative files"):
        for file_path in rep_files:
            try:
                # Get genome info
                genome = genome_lookup.get(file_path)
                if not genome:
                    logger.warning(f"No genome info found for file: {file_path}")
                    continue

                # Get destination path
                dest_path = os.path.join(out_dir, os.path.basename(file_path))

                # Copy file
                shutil.copy2(file_path, dest_path)
                logger.debug(f"Copied {file_path} to {dest_path}")

            except Exception as e:
                logger.error(f"Error copying file {file_path}: {e}")


def save_results_to_tsv(
    results: List[Dict[str, Any]], genomes: List[Dict[str, Any]], output_file: str
) -> None:
    """
    Save processing results to a TSV file.

    Args:
        results: List of processing results
        genomes: List of all genome dictionaries
        output_file: Output file path
    """
    # Create lookup for representative status
    rep_lookup = {}
    for result in results:
        if result.get("status") in ("success", "singleton"):
            for rep_file in result.get("representatives", []):
                rep_lookup[rep_file] = True

    # Create output rows
    rows = []
    for genome in genomes:
        row = {
            "taxonomy": genome["taxonomy"],
            "accession": genome["accession"],
            "file": genome["file"],
            "representative": rep_lookup.get(genome["file"], False),
        }
        rows.append(row)

    # Write to file
    with open(output_file, "w") as f:
        # Write header
        f.write("taxonomy\taccession\tfile\trepresentative\n")

        # Write rows
        for row in rows:
            f.write(
                f"{row['taxonomy']}\t{row['accession']}\t{row['file']}\t{1 if row['representative'] else 0}\n"
            )

    logger.info(f"Results saved to {output_file}")


def main():
    """Main entry point."""
    # Parse command line arguments
    args = parse_args()

    # Create configuration
    config = Config.from_args(args)

    # Initial logging setup (will be updated later with file)
    from derep_genomes.utils.logging import (
        setup_logging,
        is_debug_enabled,
        log_memory_usage,
    )

    log_level = logging.DEBUG if config.debug else logging.INFO
    setup_logging(log_level=log_level)

    logger.info(f"Starting genome dereplication with {config.similarity_method}")
    logger.info(str(config))

    # Setup run directories
    dirs = setup_run_directories(config)

    # Save configuration
    config_file = dirs["results"] / "config.yaml"
    config.save(str(config_file))
    logger.info(f"Configuration saved to {config_file}")

    # Initialize database
    db_path = dirs["results"] / f"derep-genomes_{config.timestamp}.db"
    config.db_path = str(db_path)  # Update config with actual path

    with DuckDBManager(config.db_path) as db:
        # Check or create schema
        if not check_schema(db):
            logger.info("Creating database schema")
            create_schema(db)

        # Create database operations
        db_ops = DBOperations(db)

        # Load genome data
        logger.info(f"Loading genome data from {config.tax_file}")
        genomes = load_genome_data(config.tax_file)
        logger.info(f"Loaded {len(genomes)} genomes")

        # Filter by selected taxa if provided
        if config.selected_taxa:
            logger.info(f"Filtering by selected taxa from {config.selected_taxa}")
            selected_taxa = load_taxa_from_file(config.selected_taxa)
            genomes = filter_by_taxa(genomes, selected_taxa)
            logger.info(
                f"Filtered to {len(genomes)} genomes in {len(selected_taxa)} taxa"
            )

        # Check if files exist
        logger.info("Checking if genome files exist")
        genomes = check_files_exist(genomes)
        logger.info(f"Found {len(genomes)} existing genome files")

        # Group by taxonomy
        taxonomy_groups = group_by_taxonomy(genomes)
        logger.info(f"Found {len(taxonomy_groups)} unique taxa")

        # Process taxa in parallel
        logger.info(f"Processing taxa using {config.workers} workers")

        # Create task manager
        task_manager = ParallelTaskManager(
            max_workers=config.workers, show_progress=True
        )

        # Process taxa
        results = []
        with Timer("Taxa processing"):
            for taxon, taxon_genomes in taxonomy_groups.items():
                results.append(process_taxon(taxon, taxon_genomes, config, db_ops))

        # Copy representative files if requested
        if config.copy:
            copy_representative_files(results, genomes, dirs["reps"])

        # Save results to TSV
        output_file = dirs["results"] / f"derep-genomes_results_{config.timestamp}.tsv"
        save_results_to_tsv(results, genomes, str(output_file))

        # Log summary
        n_total = len(genomes)
        n_reps = sum(
            len(result.get("representatives", []))
            for result in results
            if result.get("status") in ("success", "singleton")
        )

        logger.info(
            f"Dereplication complete. {n_reps} representatives out of {n_total} genomes ({n_reps/n_total*100:.2f}%)"
        )

        # Log memory usage
        log_memory_usage("Final memory usage")
        from derep_genomes.utils.helpers import genome_length_cache

        logger.info(f"Genome length cache stats: {genome_length_cache.get_stats()}")


if __name__ == "__main__":
    main()
