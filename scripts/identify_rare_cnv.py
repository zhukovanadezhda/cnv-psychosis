"""
A script that identifies rare CNVs from a processed VCF file.

Usage:
======
    python identify_rare_cnv.py --input <your_cnv_file> --db <db_file> --output <output_file>

Arguments:
==========
    your_cnv_file: str
        A path to the CNV file.
    db_file: str
        A path to the database with CNVs.
    output_file: str
        A path to the output CSV file.
Returns:
========
    A CSV file with rare CNVs only.
"""

__authors__ = "Nadezhda Zhukova"
__contact__ = "nadezhda.zhukova@inserm.fr"
__copyright__ = "MIT"
__date__ = "2024-05-01"
__version__ = "1.0.0"

import os
import pandas as pd
import logging
import argparse
from tqdm.auto import tqdm
tqdm.pandas()

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def get_arguments():
    """Get script arguments from command line.

    Returns
    -------
    tuple
        Filenames for CNV input, CNV database, and output file.
    """
    parser = argparse.ArgumentParser(description='Identify rare CNVs from a CNV file.')
    parser.add_argument('--input', required=True, help='Path to the CNV file.')
    parser.add_argument('--db', required=True, help='Path to the CNV database file.')
    parser.add_argument('--output', required=True, help='Path to the output CSV file.')
    args = parser.parse_args()
    return args.input, args.db, args.output

def define_rare_cnv(chromosome, start, end, cnv_type, df_all_variants, overlap_threshold=50):
    """Define rare CNVs.

    Args:
        chromosome (str): Chromosome number. e.g. "chr1".
        start (int): Start position of the CNV.
        end (int): End position of the CNV.
        cnv_type (str): Type of the CNV. e.g. "DEL".
        df_all_variants (pd.DataFrame): A DataFrame with all CNVs.
        overlap_threshold (int): The percentage of overlap between two CNVs.

    Returns:
        bool: True if the CNV is rare, False otherwise.
    """
    # Get all non-rare CNVs on the same chromosome and of the same type
    non_rare_cnvs = df_all_variants[df_all_variants["non_rare"] & 
                                    (df_all_variants["chrom"] == chromosome) & 
                                    (df_all_variants["varType"] == cnv_type)]

    if non_rare_cnvs.empty:
        logging.debug(f"No non-rare CNVs found for {chromosome}:{start}-{end} of type {cnv_type}")
        return True
    
    # Exclude the CNVs that don't overlap for sure
    non_rare_cnvs = non_rare_cnvs[(non_rare_cnvs["chromEnd"] > start) & 
                                  (non_rare_cnvs["chromStart"] < end)]
    
    if non_rare_cnvs.empty:
        logging.debug(f"No overlapping non-rare CNVs found for {chromosome}:{start}-{end} of type {cnv_type}")
        return True
    
    # Check for overlap with any non-rare CNV
    for _, row in non_rare_cnvs.iterrows():
        # Calculate the overlap between the two CNVs
        overlap_start = max(start, row["chromStart"])
        overlap_end = min(end, row["chromEnd"])
        overlap_length = max(0, overlap_end - overlap_start)
        overlap_percentage = (overlap_length / (end - start)) * 100
        
        # If the overlap is more than the threshold, the CNV is not rare
        if overlap_percentage > overlap_threshold:
            logging.debug(f"CNV {chromosome}:{start}-{end} of type {cnv_type} overlaps with non-rare CNV {row['chrom']}:{row['chromStart']}-{row['chromEnd']}")
            return False
    
    logging.debug(f"CNV {chromosome}:{start}-{end} of type {cnv_type} is rare")
    return True

if __name__ == "__main__":
    # Get script arguments
    CNV_FILENAME, DB_FILENAME, OUTPUT_FILENAME = get_arguments()
    
    logging.info(f"Reading CNV database from {DB_FILENAME}")
    # Define the column names for the CNV database
    db_columns = [
        'bin', 'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 
        'varType', 'reference', 'pubMedId', 'method', 'platform', 'mergedVariants', 'supportingVariants', 
        'sampleSize', 'observedGains', 'observedLosses', 'cohortDescription', 'genes', 'samples'
    ]
    # Read the CNV database with column names
    df_all_variants = pd.read_csv(DB_FILENAME, sep='\t', header=None, names=db_columns)
    
    logging.info("Filtering and processing CNV database")
    # Filter and process the CNV database
    columns_to_keep = ["chrom", "chromStart", "chromEnd", "varType", "sampleSize", "observedGains", "observedLosses", "genes"]
    df_all_variants = df_all_variants[columns_to_keep]
    
    cnv_types = ['gain+loss', 'loss', 'duplication', 'deletion', 'gain']
    df_all_variants = df_all_variants[df_all_variants["varType"].isin(cnv_types)]
    
    # Split 'gain+loss' into 'gain' and 'loss'
    df_gain = df_all_variants[df_all_variants["varType"] == "gain+loss"].copy()
    df_gain["varType"] = "gain"
    df_gain["observedLosses"] = 0
    
    df_loss = df_all_variants[df_all_variants["varType"] == "gain+loss"].copy()
    df_loss["varType"] = "loss"
    df_loss["observedGains"] = 0
    
    df_all_variants.drop(df_all_variants[df_all_variants["varType"] == "gain+loss"].index, inplace=True)
    df_all_variants = pd.concat([df_all_variants, df_gain, df_loss], axis=0)
    df_all_variants.sort_values(by=["chrom", "chromStart"], inplace=True)
    df_all_variants["varType"] = df_all_variants["varType"].apply(lambda x: "DUP" if x in ("gain", "duplication") else "DEL")
    df_all_variants.reset_index(drop=True, inplace=True)
    
    df_all_variants["non_rare"] = df_all_variants.apply(lambda row: (row["observedGains"] + row["observedLosses"]) > 0.01 * row["sampleSize"], axis=1)
    
    logging.info(f"Reading CNV file from {CNV_FILENAME}")
    # Read the CNV file
    df_annotated_cnv = pd.read_csv(CNV_FILENAME)
    
    logging.info("Identifying rare CNVs")
    # Define rare CNVs in the CNV file with progress bar
    df_annotated_cnv["is_rare"] = df_annotated_cnv.progress_apply(
        lambda row: define_rare_cnv(
            row["Chromosome"], row["Start"], row["End"], row["Type"], df_all_variants
        ), axis=1
    )
    
    logging.info(f"Saving results to {OUTPUT_FILENAME}")
    # Save the result to a new file
    df_annotated_cnv.to_csv(OUTPUT_FILENAME, index=False)
    logging.info(f"Rare CNVs identified and saved to {OUTPUT_FILENAME}")
