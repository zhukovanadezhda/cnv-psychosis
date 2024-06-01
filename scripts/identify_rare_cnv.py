"""

A script that identifies rare CNVs from a VCF file.

Usage:
======
    python3 identify_rare_cnv.py <your_cnv_file> <db_file>

Arguments:
==========
    your_cnv_file: str
        A path to the CNV file.
    db_file: str
        A path to the database with CNVs.
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
    if len(sys.argv) != 3:
        sys.exit("Wrong number of arguments. \n Usage: "
                 "python3 identify_rare_cnv.py <your_cnv_file> <db_file>")
    for filename in sys.argv[1:]:
        if not os.path.exists(filename):
            sys.exit(f"File {filename} does not exist. Check the path.")
    return sys.argv[1], sys.argv[2]


def define_rare_cnv(chromosome, start, end, cnv_type,
                    df_all_variants, overlap_threshold=50):
    """Define rare CNVs.

    Rare CNVs are defined as CNVs that do not overlap with any non-rare CNVs.
    Non-rare CNVs are CNVs that are observed in more than 1% of the samples.

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
    non_rare_cnvs = df_all_variants[df_all_variants["non_rare"]
                                    & (df_all_variants["chrom"] == chromosome)
                                    & (df_all_variants["varType"] == cnv_type)]
    # Exclude the CNVs that don't overlap for sure
    # (i.e. the CNVs that are completely before or after the CNV of interest)
    non_rare_cnvs = non_rare_cnvs[(non_rare_cnvs["chromEnd"] < start)
                                  | (non_rare_cnvs["chromStart"] < end)]
    # Check for overlap with any non-rare CNV
    for _, row in non_rare_cnvs.iterrows():
        # Calculate the overlap between the two CNVs
        overlap_start = max(start, row["chromStart"])
        overlap_end = min(end, row["chromEnd"])
        overlap_length = max(0, overlap_end - overlap_start)
        overlap_percentage = (overlap_length / (end - start)) * 100
        # If the overlap is more than the threshold, the CNV is not rare
        if overlap_percentage > overlap_threshold:
            return False
    return True


if __name__ == "__main__":
    # Get script arguments
    CNV_FILENAME, DB_FILENAME = get_arguments()
    # Read the CNV database
    DB_CNV = pd.read_csv(DB_FILENAME)
    # Read the CNV file
    DF_CNV = pd.read_csv(CNV_FILENAME)
    # Define non-rare CNVs in the database
    DB_CNV["non_rare"] = DB_CNV.apply(
        lambda row: (
            row["observedGains"] + row["observedLosses"]
        ) > 0.01 * row["sampleSize"], axis=1
        )
    # Define rare CNVs in the CNV file
    DB_CNV["non_rare"] = DB_CNV.apply(
        lambda row: (
            row["observedGains"] + row["observedLosses"]
            ) > 0.01 * row["sampleSize"], axis=1
        )
    print(DF_CNV["is_rare"].value_counts())
    print(f"Your updated file is saved as "
          f"{CNV_FILENAME.split('.')[0]}_rare.csv")
    # Save the result to a new file
    DF_CNV.to_csv(f"{CNV_FILENAME.split('.')[0]}_rare.csv", index=False)
