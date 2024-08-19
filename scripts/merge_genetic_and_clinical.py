"""
Script to merge clinical data with CNV statistics.

This script merges clinical data from a metadata CSV file with the CNV statistics
generated for each individual based on their ID.

Usage:
======
    python merge_genetic_and_clinical.py <cnv_statistics_csv> <output_csv_file>

Arguments:
==========
    cnv_statistics_csv : str
        Path to the input CSV file containing per-individual CNV statistics.
    output_csv_file : str
        Path to save the merged CSV file with clinical data.
"""

import sys
import logging
import pandas as pd


def setup_logging() -> None:
    """Set up the logging configuration."""
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def merge_clinical_and_genetical_data(cnv_statistics_csv: str, 
                                      metadata_csv: str, 
                                      output_csv_file: str) -> None:
    """
    Merge clinical data from metadata CSV with CNV statistics.

    Parameters:
    ===========
    cnv_statistics_csv : str
        Path to the input CSV file with per-individual CNV statistics.
    metadata_csv : str
        Path to the metadata CSV file with clinical data.
    output_csv_file : str
        Path to save the merged CSV file.
    """
    try:
        logging.info("Loading CNV statistics data from %s", cnv_statistics_csv)
        df_cnv_stats = pd.read_csv(cnv_statistics_csv)

        logging.info("Loading clinical metadata from %s", metadata_csv)
        df_metadata = pd.read_csv(metadata_csv)

        logging.info("Merging CNV statistics with clinical metadata on 'ID'")
        df_merged = pd.merge(df_cnv_stats, df_metadata, on="ID", how="left")

        logging.info("Saving merged data to %s", output_csv_file)
        df_merged.to_csv(output_csv_file, index=False)

        logging.info("Merging completed successfully")

    except Exception as e:
        logging.error("Error during merging: %s", e)
        sys.exit(1)

if __name__ == "__main__":
    setup_logging()

    if len(sys.argv) != 3:
        print("Usage: python merge_genetical_and_clinical.py "
              "<cnv_statistics_csv> <output_csv_file>")
        sys.exit(1)

    cnv_statistics_csv = sys.argv[1]
    metadata_csv = "data/metadata.csv"
    output_csv_file = sys.argv[2]

    merge_clinical_and_genetical_data(cnv_statistics_csv,
                                      metadata_csv, 
                                      output_csv_file)
