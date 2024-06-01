"""

A script that converts gene names to Ensembl IDs and vice versa.

Usage:
======
    python3 convert_gene_id.py <ensemble_to_name=True> <your_id_list>

Arguments:
==========
    ensemble_to_name: str (True or False)
        If True, the script converts Ensembl IDs to gene names.
        If False, the script converts gene names to Ensembl IDs.
    your_id_list: str (comma-separated, no spaces)
        A list of gene names or Ensembl IDs to convert.
Returns:
========
    Converted list of gene names or Ensembl IDs.
"""

__authors__ = "Nadezhda Zhukova"
__contact__ = "nadezhda.zhukova@inserm.fr"
__copyright__ = "MIT"
__date__ = "2024-05-01"
__version__ = "1.0.0"

import json
import os
import pandas as pd
import sys


def get_arguments():
    """Get script arguments from command line.

    Returns
    -------
    str
        The type of conversion (Ensembl ID to gene name or vice versa).
    str
        A csv file with gene names or Ensembl IDs to convert.
    """
    # Check the number of arguments
    if len(sys.argv) != 3:
        sys.exit("Wrong number of arguments. \n Usage: "
                 "python3 convert_gene_id.py <True/False> <your_csv_path>")
    # Check if the file exists
    if not os.path.exists(sys.argv[2]):
        sys.exit(f"File {sys.argv[2]} does not exist.")
    # Check the type of conversion
    if sys.argv[1] not in ['True', 'False']:
        sys.exit("The second argument should be either True "
                 "(for Ensembl ID to gene name converter) or "
                 "False (for gene name to Ensembl ID converter).")
    return sys.argv[1], sys.argv[2]


def convert_gene_id(csv_path, json_path):
    """Convert gene names to Ensembl IDs and vice versa.

    Args:
        csv_path (str): Path to the CSV file with gene names.
        json_path (str): Path to the JSON file with gene names and Ensembl IDs.

    Returns:
        pd.DataFrame: DataFrame with gene names and Ensembl IDs.
    """
    with open(json_path, 'r', encoding='utf-8') as f_in:
        gene_ids = json.load(f_in)
    gene_ids_first = {k: v.split(", ")[0] for k, v in gene_ids.items()}
    df = pd.read_csv(csv_path)
    df["Gene_name"] = df["Gene"].apply(lambda x: gene_ids_first.get(x, None))
    missing_genes = [gene for gene, ensembl_id
                     in zip(df["Gene"].to_list(), 
                            df["Gene_name"].to_list()) 
                     if ensembl_id is None]
    if missing_genes:
        print(f"No data found for genes: {', '.join(missing_genes)}")
    return df


if __name__ == "__main__":
    # Get script arguments
    CONVERTER_TYPE, CSV_PATH = get_arguments()
    # Choose the JSON file according to the type of conversion
    if CONVERTER_TYPE == 'True':
        JSON_PATH = "data/ensemble_to_name.json"
    else:
        JSON_PATH = "data/name_to_ensemble.json"
    # Check if the JSON file exists
    if not os.path.exists(JSON_PATH):
        sys.exit(f"File {JSON_PATH} does not exist. "
                 "JSON file should be in the same directory with this script.")
    # Perform the conversion
    converted_csv = convert_gene_id(CSV_PATH, JSON_PATH)
    # Save the new file
    new_name = CSV_PATH.split(".")[0] + "_converted.csv"
    converted_csv.to_csv(new_name, index=False)
