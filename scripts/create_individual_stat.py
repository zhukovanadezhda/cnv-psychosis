"""
Script to generate CNV statistics per individual.

This script processes a CSV file containing annotated CNVs and produces
a summary of CNV statistics for each individual, including counts, sums,
and classifications.

Usage:
======
    python create_stat_per_individual.py <input_csv_file> <output_csv_file>

Arguments:
==========
    input_csv_file : str
        Path to the input CSV file containing annotated CNVs.
    output_csv_file : str
        Path to save the output CSV file with per-individual statistics.
"""

import sys
import logging
import pandas as pd
from tqdm.auto import tqdm

def setup_logging() -> None:
    """Set up the logging configuration."""
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def compute_statistics_per_individual(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute statistics per individual from the annotated CNV DataFrame.

    Parameters:
    ===========
    df : pd.DataFrame
        DataFrame containing annotated CNV data.

    Returns:
    ========
    pd.DataFrame
        DataFrame containing per-individual CNV statistics.
    """
    logging.info("Computing statistics per individual")

    # Fill missing values for gene columns
    df["All protein coding genes"] = df["All protein coding genes"].fillna("")
    df["Brain_genes"] = df["Brain_genes"].fillna("")
    df["Known or predicted dosage-sensitive genes"] = df["Known or predicted dosage-sensitive genes"].fillna("")
    df["is_long"] = (df["End"] - df["Start"]) > 1_000_000


    # Calculate additional flags
    df["is_long_brain"] = df["Is_brain"] & df["is_long"]

    # Group and aggregate data
    aggregated = df.groupby(["ID", "Type"]).agg(
        CNV_length=("End", lambda x: (x - df.loc[x.index, "Start"]).sum()),
        NB_CNV=("ID", "size"),
        SUM_quality=("Quality", "sum"),
        Copy_number=("Copy_number", "sum"),        
        NB_rare=("is_rare", "sum"),
        NB_long=("is_long", "sum"),
        NB_brain=("Is_brain", "sum"),
        NB_long_rare=("is_long_brain", "sum"),
        Rare_length=("End", lambda x: (x[df["is_rare"]] - df.loc[x[df["is_rare"]].index, "Start"]).sum()),
        Protein_coding_genes=("All protein coding genes", lambda x: ', '.join(set(x))),
        Brain_genes=("Brain_genes", lambda x: ', '.join(set(x))),
        Dosage_sensitive_genes=("Known or predicted dosage-sensitive genes", lambda x: ', '.join(set(x))),
        Likely_benign=("Classification", lambda x: (x == "Likely benign").sum()),
        Likely_pathogenic=("Classification", lambda x: (x == "Likely pathogenic").sum()),
        Pathogenic=("Classification", lambda x: (x == "Pathogenic").sum()),
        Uncertain_significance=("Classification", lambda x: (x == "Uncertain significance").sum())
        
    ).reset_index()

    # Calculate gene counts and rare CNV gene counts
    aggregated["NB_genes"] = aggregated["Protein_coding_genes"].apply(lambda row: len(row.split(", ")) if row else 0)
    aggregated["NB_brain_genes"] = aggregated["Brain_genes"].apply(lambda row: len(row.split(", ")) if row else 0)

    df_rare = df[df["is_rare"]].groupby(["ID", "Type"]).agg(
        Protein_coding_genes_in_rare_CNV=("All protein coding genes", lambda x: ', '.join(set(x))),
        Brain_genes_in_rare_CNV=("Brain_genes", lambda x: ', '.join(set(x))),
    ).reset_index()

    df_rare["NB_genes_rare"] = df_rare["Protein_coding_genes_in_rare_CNV"].apply(lambda row: len(row.split(", ")) if row else 0)
    df_rare["NB_brain_genes_rare"] = df_rare["Brain_genes_in_rare_CNV"].apply(lambda row: len(row.split(", ")) if row else 0)

    # Merge rare CNV data back into the main DataFrame
    aggregated = aggregated.merge(df_rare, on=["ID", "Type"], how="left")

    # Fill missing values in the merged DataFrame
    aggregated.fillna({
        "Rare_length": 0,
        "NB_genes_rare": 0,
        "NB_brain_genes_rare": 0,
        "Protein_coding_genes_in_rare_CNV": "",
        "Brain_genes_in_rare_CNV": ""
    }, inplace=True)

    # Map CNV type to a numeric value
    aggregated["Type_numeric"] = aggregated["Type"].map({"DEL": 1, "DUP": 2})

    return aggregated

def main(input_csv_file: str, output_csv_file: str) -> None:
    """
    Main function to generate per-individual statistics of CNVs.

    Parameters:
    ===========
    input_csv_file : str
        Path to the input CSV file with annotated CNVs.
    output_csv_file : str
        Path to save the output CSV file with per-individual statistics.
    """
    try:
        logging.info("Loading annotated CNV data from %s", input_csv_file)
        df_annotated_cnv = pd.read_csv(input_csv_file)

        logging.info("Computing statistics for each individual")
        tqdm.pandas(desc="Processing individuals")
        df_statistics = compute_statistics_per_individual(df_annotated_cnv)

        logging.info("Saving per-individual statistics to %s", output_csv_file)
        df_statistics.to_csv(output_csv_file, index=False)

        logging.info("Processing completed successfully")

    except Exception as e:
        logging.error("Error during processing: %s", e)
        sys.exit(1)

if __name__ == "__main__":
    setup_logging()
    
    if len(sys.argv) != 3:
        print("Usage: python create_stat_per_individual.py "
              "<input_csv_file> <output_csv_file>")
        sys.exit(1)

    input_csv_file = sys.argv[1]
    output_csv_file = sys.argv[2]

    main(input_csv_file, output_csv_file)
