"""
Filter out pathogenic CNVs from a dataset of CNV annotations.

This script takes two input CSV files:
1. A file containing all CNV annotations.
2. A file containing pathogenic CNVs.

It outputs a filtered CSV file with the pathogenic CNVs removed, 
based on the common index columns 'ID' and 'Cytoband'.

Usage:
    python remove_pathogenic_cnv.py <cnv_file> <patho_file> <output_file>

Arguments:
    cnv_file    : Path to the input CNV annotations file.
    patho_file  : Path to the input pathogenic CNV file.
    output_file : Path to the output filtered CNV annotations file.

Example:
    python remove_pathogenic_cnv.py \ 
    results/merged_all_cnv_annotations.csv \
    results/pathogenic_cnv_with_marshal.csv \
    results/uncertain_merged_all_cnv_annotations_filtered.csv
"""

import pandas as pd
import argparse
import logging


def setup_logging():
    """Configure logging to display timestamps and log levels."""
    logging.basicConfig(
        level=logging.INFO, 
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def load_dataframe(filepath, index_columns):
    """
    Load a CSV file into a DataFrame and set the specified columns as the index.
    
    Parameters:
    filepath (str): Path to the CSV file.
    index_columns (list): List of columns to set as the index.
    
    Returns:
    pd.DataFrame: The loaded DataFrame with the index set.
    """
    logging.info(f"Loading data from {filepath}")
    df = pd.read_csv(filepath)
    df.set_index(index_columns, inplace=True)
    return df


def filter_pathogenic_cnv(df_all, df_pathogenic):
    """
    Remove rows corresponding to pathogenic CNVs from the full dataset.
    
    Parameters:
    df_all (pd.DataFrame): DataFrame containing all CNV annotations.
    df_pathogenic (pd.DataFrame): DataFrame containing pathogenic CNV annotations.
    
    Returns:
    pd.DataFrame: The filtered DataFrame with pathogenic CNVs removed.
    """
    logging.info("Filtering pathogenic CNVs from the dataset")
    indices_to_remove = df_pathogenic.index
    filtered_df = df_all.loc[df_all.index.difference(indices_to_remove)]
    return filtered_df


def save_dataframe(df, output_filepath):
    """
    Save the DataFrame to a CSV file.
    
    Parameters:
    df (pd.DataFrame): The DataFrame to be saved.
    output_filepath (str): Path to the output CSV file.
    """
    logging.info(f"Saving filtered data to {output_filepath}")
    df.to_csv(output_filepath)


def main(cnv_file, patho_file, output_file):
    """
    Main function to execute the CNV filtering process.
    
    Parameters:
    cnv_file (str): Path to the input CNV annotations file.
    patho_file (str): Path to the input pathogenic CNV file.
    output_file (str): Path to the output filtered CNV annotations file.
    """
    # Initialize logging
    setup_logging()

    # Load the full CNV dataset and the pathogenic CNV dataset
    df_all = load_dataframe(cnv_file, ['ID', 'Cytoband'])
    df_pathogenic = load_dataframe(patho_file, ['ID', 'Cytoband'])

    # Filter out pathogenic CNVs
    filtered_df = filter_pathogenic_cnv(df_all, df_pathogenic)

    # Save the filtered dataset
    save_dataframe(filtered_df, output_file)


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Remove pathogenic CNVs from a dataset."
    )
    parser.add_argument('cnv_file', 
                        help="Path to the input CNV annotations file.")
    parser.add_argument('patho_file', 
                        help="Path to the input pathogenic CNV file.")
    parser.add_argument('output_file', 
                        help="Path to the output filtered CNV annotations file.")
    
    # Get the arguments from the command line
    args = parser.parse_args()

    # Execute the main function with parsed arguments
    main(args.cnv_file, args.patho_file, args.output_file)
