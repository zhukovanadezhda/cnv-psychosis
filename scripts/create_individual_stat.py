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


def aggregate_group(group):

    result = {
        "Avg_quality": group["Quality"].mean(),
        "Avg_copy_number": group["Copy_number"].mean()
    }

    flags = [
        "Is_rare",
        "Is_long_rare",
        "Is_brain_rare",
        "Is_only_brain_rare",
        "Is_long_brain_rare",
        "Is_long_only_brain_rare"
    ]
    
    # Count occurrences and lengths
    for flag in flags:
        # Remove 'Is_' prefix
        suffix = flag.lower()[3:]
        mask = group[flag]
        result[f"NB_{suffix}"] = mask.sum()
        result[f"{suffix.capitalize()}_length"] = (
            group.loc[mask, "End"] - group.loc[mask, "Start"]
            ).sum()

    # Genes (unique concatenation)
    gene_fields = {
        "Protein_coding_genes": "All protein coding genes",
        "Brain_genes": "Brain_genes",
        "Only_brain_genes": "Only_brain_genes",
        "Dosage_sensitive_genes": "Known or predicted dosage-sensitive genes"
    }
    for key, col in gene_fields.items():
        # Concatenate unique gene names
        result[key] = ', '.join(sorted(set(group[col].dropna())))
        # Count number of unique genes in rare CNVs
        rare_mask = group["Is_rare"] & group[col].notna()
        result[f"NB_{key.lower()}_rare"] = len(set(group.loc[rare_mask, col]))

    return pd.Series(result)


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

    # Calculate additional flags
    # At this point, we already have 'Is_brain', 'Is_only_brain', 'Is_rare'
    df["Is_long"] = (df["End"] - df["Start"]) > 1_000_000
    df["Is_long_brain"] = df["Is_brain"] & df["Is_long"]
    df["Is_long_only_brain"] = df["Is_only_brain"] & df["Is_long"]
    df["Is_long_rare"] = df["Is_rare"] & df["Is_long"]
    df["Is_brain_rare"] = df["Is_brain"] & df["Is_rare"]
    df["Is_only_brain_rare"] = df["Is_only_brain"] & df["Is_rare"]
    df["Is_long_brain_rare"] = df["Is_long_brain"] & df["Is_rare"]
    df["Is_long_only_brain_rare"] = df["Is_long_only_brain"] & df["Is_rare"]

    # Group by 'ID' and 'Type' and apply the aggregation function
    df = df.groupby(["ID", "Type"]).apply(aggregate_group).reset_index()

    # Pivot the table to create separate columns for 'del' and 'dup'
    df_pivot = df.pivot_table(
        index='ID',
        columns='Type',
        values=[col for col in df.columns if col not in ["ID", "Type"]],
        aggfunc='first'
    )

    # Flatten columns
    df_pivot.columns = ['_'.join(col).strip() for col in df_pivot.columns.values]
    df_pivot.reset_index(inplace=True)

    return df_pivot

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
