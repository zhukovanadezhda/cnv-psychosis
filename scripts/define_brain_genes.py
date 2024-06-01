"""

A script that defines brain genes from a VCF file.

Usage:
======
    python3 define_brain_genes.py <your_vcf_file> <db_brain>

Arguments:
==========
    your_vcf_file: str
        A path to the VCF file.
    db_brain: str
        A path to the database with brain genes.
Returns:
========
    A CSV file with brain SNP/CNV only.
"""

__authors__ = "Nadezhda Zhukova"
__contact__ = "nadezhda.zhukova@inserm.fr"
__copyright__ = "MIT"
__date__ = "2024-05-01"
__version__ = "1.0.0"

import os
import sys
import pandas as pd
from tqdm.auto import tqdm
tqdm.pandas()


def get_arguments():
    """Get script arguments from command line.

    Returns
    -------
    str
        Filenames.
    """
    # Check the number of arguments
    if len(sys.argv) != 3:
        sys.exit("Wrong number of arguments. \n Usage: "
                 "python3 define_brain_genes.py <your_vcf_file> <db_brain>")
    # Check the existence of files
    for filename in sys.argv[1:]:
        if not os.path.exists(filename):
            sys.exit(f"File {filename} does not exist. Check the path.")
    return sys.argv[1], sys.argv[2]


def get_vcf_names(vcf_path):
    """Get column names from the VCF file.

    Args:
        vcf_path (str): Path to the VCF file.

    Returns:
        list: A list of column names from the VCF file.
    """
    with open(vcf_path, "r", encoding="utf-8") as in_file:
        for line in in_file:
            # Get the header
            if line.startswith("#Uploaded_variation"):
                vcf_names = [x.strip() for x in line.split('\t')]
                break
    return vcf_names


def process_vcf(vcf_path, brain_genes_list, chunksize=100_000):
    """Process the VCF file and define brain genes.

    Args:
        vcf_path (str): Path to the VCF file.
        brain_genes_list (list): A list of brain genes.
        chunksize (int): The number of rows to read at a time. Default 100000.

    Returns:
        DataFrame: A DataFrame with brain genes only.
    """
    # Get column names from the VCF file
    names = get_vcf_names(vcf_path)
    # Define the type of columns (to save memory usage)
    dtype = {col: str for col in names}
    df_vcf = pd.DataFrame()
    # Read the VCF file by chunks
    chunks = pd.read_csv(vcf_path, comment='#', sep="\t", header=None,
                         names=names, dtype=dtype, chunksize=chunksize)
    # Process each chunk
    for i, chunk in enumerate(chunks):
        print(f"Processing chunk {i+1}: "
              f"from {chunk.index[0]} to {chunk.index[-1]}")
        # Define brain genes in the chunk
        chunk["Is_brain"] = chunk["Gene"].apply(lambda x: x in brain_genes_list)
        # Save brain genes to the DataFrame
        df_vcf = pd.concat([df_vcf, chunk[chunk["Is_brain"]]])
    return df_vcf


if __name__ == "__main__":
    # Get script arguments
    VCF_FILENAME, DB_FILENAME = get_arguments()
    # Read the database with brain genes
    DF_BRAIN = pd.read_csv(DB_FILENAME, sep="\t")
    # Get the list of brain genes
    LIST_BRAIN = DF_BRAIN["Ensembl"].to_list()
    print(f"{len(LIST_BRAIN)} brain genes found in your database.")
    print(f"Processing {VCF_FILENAME} file...")
    print("File will be processed by chunks. Please wait...")
    # Process the VCF file
    DF_VCF = process_vcf(VCF_FILENAME, LIST_BRAIN)
    # Save the result to a new file
    NEW_FILENAME = f"{VCF_FILENAME.split('.')[0]}_brain.csv"
    DF_VCF.to_csv(NEW_FILENAME, index=False)
    print(f"{len(DF_VCF)} brain SNP/CNV found.")
    print(f"{len(list(DF_VCF['Gene'].unique()))} brain SNP/CNV found.")
    print(f"The file with brain modifications only is saved as {NEW_FILENAME}")
