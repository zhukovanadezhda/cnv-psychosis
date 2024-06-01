"""A script creating two dictionaries mapping Ensembl IDs to gene names.

Usage:
======
    python3 create_gene_id_mapping.py

Returns:
========
    Two JSON files with dictionaries mapping Ensembl IDs to gene names.
"""

__authors__ = "Nadezhda Zhukova"
__contact__ = "nadezhda.zhukova@inserm.fr"
__copyright__ = "MIT"
__date__ = "2024-05-01"
__version__ = "1.0.0"

import json
import pandas as pd


def create_ensemble_to_name_mapping(input_file, output_file):
    """Create and save a dictionary mapping Ensemble IDs to gene names.

    Args:
        input_file (str): Path to the CSV file with gene data.
        output_file (str): Path to the output JSON file.

    Returns:
        dict: A dictionary mapping Ensemble IDs to gene names.
    """
    # Read data from CSV file into a DataFrame
    df_genes = pd.read_csv(input_file, sep="\t")

    # Fill missing values with empty strings
    df_genes['Gene name'].fillna("", inplace=True)
    df_genes['Gene Synonym'].fillna("", inplace=True)

    # Concatenate 'Gene name' and 'Gene Synonym' to create 'Gene names' column
    df_genes['Gene names'] = (
        df_genes['Gene name'] + ', ' + df_genes['Gene Synonym']
        )

    # Group by 'Gene stable ID' and concatenate 'Gene names' for each group
    df_genes = df_genes.groupby("Gene stable ID")['Gene names'].apply(
        ', '.join
        ).reset_index(name='Gene names')

    # Strip leading/trailing whitespaces and commas from 'Gene names' column
    df_genes['Gene names'] = df_genes['Gene names'].str.strip(', ')

    # Create dictionary mapping Ensemble IDs to gene names
    ensemble_to_name = {gene_id: gene_names for gene_id, gene_names in zip(
        df_genes["Gene stable ID"], df_genes["Gene names"]
        ) if gene_names}

    # Save dictionary to JSON file
    with open(output_file, 'w', encoding="utf-8") as f_out:
        json.dump(ensemble_to_name, f_out)

    return ensemble_to_name


def create_name_to_ensemble_mapping(input_file, output_file):
    """Create and save a dictionary mapping gene names to Ensemble IDs.

    Args:
        input_file (str): Path to the CSV file with gene data.
        output_file (str): Path to the output JSON file.

    Returns:
        dict: A dictionary mapping gene names to Ensemble IDs.
    """
    # Read data from CSV file into a DataFrame
    df_genes = pd.read_csv(input_file, sep="\t")

    # Select relevant columns
    df_genes = df_genes[["Gene stable ID", "Gene name", "Gene Synonym"]]

    # Fill missing values with empty strings
    df_genes.fillna('', inplace=True)

    # Concatenate 'Gene name' and 'Gene Synonym' for each group
    df_genes = df_genes.groupby(["Gene stable ID", "Gene name"]).apply(
        lambda x: ', '.join(x["Gene Synonym"])
        ).reset_index(name='Gene Synonym')
    df_genes = df_genes.groupby("Gene stable ID").apply(
        lambda x: ', '.join(x["Gene name"] + ', ' + x["Gene Synonym"])
        ).reset_index(name='Gene names')

    # Strip leading/trailing commas and whitespaces from 'Gene names' column
    df_genes["Gene names"] = df_genes["Gene names"].str.strip(", ")

    # Filter out rows where 'Gene names' is not empty
    df_genes = df_genes[df_genes["Gene names"] != ""]

    # Create dictionary mapping gene names to Ensemble IDs
    name_to_ensemble = {name: stable_id for stable_id, names in zip(
        df_genes['Gene stable ID'],
        df_genes['Gene names'].str.split(', ')
        ) for name in names}

    # Save dictionary to JSON file
    with open(output_file, 'w', encoding="utf-8") as f_out:
        json.dump(name_to_ensemble, f_out)

    return name_to_ensemble


if __name__ == "__main__":
    DB_FILE = 'data/mart_export.txt'
    ENSEMBLE_TO_NAME_FILE = 'data/ensemble_to_name.json'
    NAME_TO_ENSEMBLE_FILE = 'data/name_to_ensemble.json'

    # Create Ensemble to Name mapping
    ENSEMBLE_TO_NAME_DICT = create_ensemble_to_name_mapping(
        DB_FILE,
        ENSEMBLE_TO_NAME_FILE
        )
    print(
        "Ensemble to Name mapping created. Length:",
        len(ENSEMBLE_TO_NAME_DICT)
        )

    # Create Name to Ensemble mapping
    NAME_TO_ENSEMBLE_DICT = create_name_to_ensemble_mapping(
        DB_FILE,
        NAME_TO_ENSEMBLE_FILE
        )
    print(
        "Name to Ensemble mapping created. Length:",
        len(NAME_TO_ENSEMBLE_DICT)
        )
