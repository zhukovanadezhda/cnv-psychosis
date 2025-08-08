"""
A script that defines brain genes from a VCF file.

Usage:
======
    python3 define_brain_genes.py --input <your_vcf_file> --db <db_brain> --output <output_file>

Arguments:
==========
    your_vcf_file: str
        A path to the VCF file.
    db_brain: str
        A path to the database with brain genes.
    output_file: str
        A path to the output CSV file.
Returns:
========
    A CSV file with brain SNP/CNV only.
"""

__authors__ = "Nadezhda Zhukova"
__contact__ = "nadezhda.zhukova@inserm.fr"
__copyright__ = "MIT"
__date__ = "2024-05-01"
__version__ = "1.0.0"

# Import necessary libraries
import argparse
import pandas as pd
import logging
from tqdm.auto import tqdm
tqdm.pandas()

# Set up logging
logging.basicConfig(
    level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s'
)


def get_arguments():
    """Get script arguments from command line.

    Returns
    -------
    tuple
        Filenames for VCF input, brain gene database, and output file.
    """
    parser = argparse.ArgumentParser(
        description='Process VCF file to find brain genes.'
        )
    parser.add_argument(
        '--input',
        required=True,
        help='Path to the VCF file.'
        )
    parser.add_argument(
        '--db-brain',
        dest='db_brain',
        required=True,
        help='Path to the brain gene database file.'
        )
    parser.add_argument(
        '--db-brain-only',
        dest='db_brain_only',
        required=True,
        help='Path to the only brain gene database file.'
        )
    parser.add_argument(
        '--output',
        required=True,
        help='Path to the output CSV file.'
        )
    args = parser.parse_args()
    return args.input, args.db_brain, args.db_brain_only, args.output


def get_vcf_names(vcf_path):
    """Get column names from the VCF file.

    Args:
        vcf_path (str): Path to the VCF file.

    Returns:
        list: A list of column names from the VCF file.
    """
    vcf_df = pd.read_csv(vcf_path, sep='\t', comment='#', nrows=1)
    return vcf_df.columns.to_list()


def process_vcf(vcf_path, brain_genes_df, only_brain_genes, chunksize=100_000):
    """Process the VCF file and define brain genes.

    Args:
        vcf_path (str): Path to the VCF file.
        brain_genes_df (DataFrame): DataFrame containing brain gene information.
        chunksize (int): The number of rows to read at a time. Default 100000.

    Returns:
        DataFrame: A DataFrame with brain genes only.
    """
    
    df_vcf = pd.DataFrame()
    
    # Read the VCF file by chunks
    chunks = pd.read_csv(vcf_path, chunksize=chunksize)
    
    # Create a dictionary for quick lookup of brain genes
    brain_genes_dict = {
        row['Gene.name']: row for idx, row in brain_genes_df.iterrows()
        }
    only_brain_genes_dict = {
        row['Gene.name']: row for idx, row in only_brain_genes.iterrows()
        }
    
    gene_column='All protein coding genes'

    # Process each chunk
    for i, chunk in enumerate(chunks):
        logging.info(f"Processing chunk {i+1}: "
                     f"from {chunk.index[0]} to {chunk.index[-1]}")
        # Define brain genes in the chunk
        chunk['Is_brain'] = chunk[gene_column].apply(
            lambda x: any(
                gene in brain_genes_dict for gene in str(x).split(", ")
            ) 
            if pd.notna(x) else False
        )
        # Define only brain genes in the chunk
        chunk['Is_only_brain'] = chunk[gene_column].apply(
            lambda x: any(
                gene in only_brain_genes_dict for gene in str(x).split(", ")
            ) 
            if pd.notna(x) else False
        )
        # Keep tracks of brain genes
        chunk['Brain_genes'] = chunk[gene_column].apply(
            lambda x: ", ".join(
                gene for gene in str(x).split(", ") if gene in brain_genes_dict
                ) if pd.notna(x) else ""
            )
        # Keep tracks of only brain genes
        chunk['Only_brain_genes'] = chunk[gene_column].apply(
            lambda x: ", ".join(
                gene for gene in str(x).split(", ") if gene in only_brain_genes_dict
                ) if pd.notna(x) else ""
            )
        # Save the DataFrame
        df_vcf = pd.concat([df_vcf, chunk])
    return df_vcf


if __name__ == "__main__":
    # Get script arguments
    VCF_FILENAME, DB_BRAIN, DB_BRAIN_ONLY, OUTPUT_FILENAME = get_arguments()
    
    # Read the database with brain genes
    DF_BRAIN = pd.read_csv(DB_BRAIN, sep="\t")
    DF_BRAIN_ONLY = pd.read_csv(DB_BRAIN_ONLY, sep="\t")
    
    logging.info(f"{len(DF_BRAIN)} brain genes found in your database.")
    logging.info(f"{len(DF_BRAIN_ONLY)} only brain genes found in your database.")
    logging.info(f"Processing {VCF_FILENAME} file...")
    logging.info("File will be processed by chunks. Please wait...")
    
    # Process the VCF file
    DF_VCF = process_vcf(VCF_FILENAME, DF_BRAIN, DF_BRAIN_ONLY)
    
    # Save the result to a new file
    DF_VCF.to_csv(OUTPUT_FILENAME, index=False)
    logging.info(f"{len(DF_VCF)} brain SNP/CNV found.")
    logging.info(f"{len(list(DF_VCF['Brain_genes'].unique()))} brain genes found.")
    logging.info(f"{len(list(DF_VCF['Only_brain_genes'].unique()))} only brain genes found.")
    logging.info(f"The file with brain modifications only is saved as {OUTPUT_FILENAME}")
