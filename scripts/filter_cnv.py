"""
A script to filter CNVs in a CSV file.

It filters out CNVs based on various criteria provided by the user.

Usage:
======
    python filter_cnv.py --csv_file input.csv --output_file output.csv 
                        [--verbose] [--cnvLength] [--cnvQual] 
                        [--cnvBinSupportRatio] [--cnvCopyRatio] [--Chromosome] 
                        [--Classification] [--is_rare] [--is_long] [--is_brain]

"""

__authors__ = "Nadezhda Zhukova"
__contact__ = "nadezhda.zhukova@inserm.fr"
__copyright__ = "MIT"
__date__ = "2024-07-29"
__version__ = "1.0.0"

import argparse
import pandas as pd
import logging


def setup_logging(verbose: bool) -> None:
    """Configure the logging settings."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def read_csv(csv_file: str) -> pd.DataFrame:
    """Read the CSV file into a DataFrame."""
    return pd.read_csv(csv_file)


def apply_filters(df: pd.DataFrame, filters: dict) -> pd.DataFrame:
    """Apply the specified filters to the DataFrame."""
    initial_count = len(df)
    logging.info(f"Initial CNV number {initial_count}.")
    if filters['cnvLength']:
        df = df[~df['Filter'].str.contains('cnvLength')]
        logging.info(f"Filtered out {initial_count - len(df)} "
                     f"CNVs due to cnvLength.")
        initial_count = len(df)
    if filters['cnvQual']:
        df = df[~df['Filter'].str.contains('cnvQual')]
        logging.info(f"Filtered out {initial_count - len(df)} "
                     f"CNVs due to cnvQual.")
        initial_count = len(df)
    if filters['cnvBinSupportRatio']:
        df = df[~df['Filter'].str.contains('cnvBinSupportRatio')]
        logging.info(f"Filtered out {initial_count - len(df)} "
                     f"CNVs due to cnvBinSupportRatio.")
        initial_count = len(df)
    if filters['cnvCopyRatio']:
        df = df[~df['Filter'].str.contains('cnvCopyRatio')]
        logging.info(f"Filtered out {initial_count - len(df)} "
                     f"CNVs due to cnvCopyRatio.")
        initial_count = len(df)
    if filters['Chromosome']:
        df = df[df['Chromosome'] != 'chrY']
        logging.info(f"Filtered out {initial_count - len(df)} "
                     f"CNVs due to Chromosome being chrY.")
        initial_count = len(df)
    if filters['Classification']:
        df = df[df['Classification'] != 'Benign']
        logging.info(f"Filtered out {initial_count - len(df)} "
                     f"CNVs due to Classification being Benign.")
        initial_count = len(df)
    if filters['is_rare'] is not None:
        df = df[df['is_rare'] == filters['is_rare']]
        logging.info(f"Filtered out {initial_count - len(df)} "
                     f"CNVs due to is_rare being {filters['is_rare']}.")
        initial_count = len(df)
    if filters['is_long'] is not None:
        df = df[df['is_long'] == filters['is_long']]
        logging.info(f"Filtered out {initial_count - len(df)} "
                     f"CNVs due to is_long being {filters['is_long']}.")
        initial_count = len(df)
    if filters['is_brain'] is not None:
        df = df[df['is_brain'] == filters['is_brain']]
        logging.info(f"Filtered out {initial_count - len(df)} "
                     f"CNVs due to is_brain being {filters['is_brain']}.")
        initial_count = len(df)
        
    logging.info(f"Total CNV number {initial_count}.")
    return df


def save_filtered_csv(df: pd.DataFrame, output_file: str) -> None:
    """Save the filtered DataFrame to a CSV file."""
    df.to_csv(output_file, index=False)


def main(csv_file: str, filters: dict, output_file: str, verbose: bool) -> None:
    """Main function to read, filter, and save the CSV file."""
    setup_logging(verbose)
    
    logging.info("Reading CSV file.")
    df = read_csv(csv_file)
    
    logging.info("Applying filters.")
    filtered_df = apply_filters(df, filters)
    
    logging.info("Saving filtered CSV file.")
    save_filtered_csv(filtered_df, output_file)
    
    logging.info("Filtering completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter CNVs in a CSV file')
    parser.add_argument('--csv_file', type=str, required=True, 
                        help='Path to the input CSV file.')
    parser.add_argument('--output_file', type=str, required=True, 
                        help='Path to the output CSV file.')
    parser.add_argument('--verbose', action='store_true', 
                        help='Enable verbose mode for logging.')
    parser.add_argument('--cnvLength', action='store_true', 
                        help='Filter out CNVs with cnvLength < 1kbp.')
    parser.add_argument('--cnvQual', action='store_true', 
                        help='Filter out CNVs with cnvQual < 10.')
    parser.add_argument('--cnvBinSupportRatio', action='store_true', 
                        help='Filter out CNVs with cnvBinSupportRatio < 0.2.')
    parser.add_argument('--cnvCopyRatio', action='store_true', 
                        help='Filter out CNVs with cnvCopyRatio out of 1+-0.2')
    parser.add_argument('--Chromosome', action='store_true',
                        help='Filter out CNVs on the Y chromosome.')
    parser.add_argument('--Classification', action='store_true',
                        help='Filter out CNVs classified as Benign.')
    parser.add_argument('--is_rare', type=bool,
                        help='Filter for rare CNVs (True/False).')
    parser.add_argument('--is_long', type=bool,
                        help='Filter for long CNVs (True/False).')
    parser.add_argument('--is_brain', type=bool,
                        help='Filter for brain CNVs (True/False).')

    args = parser.parse_args()

    filters = {
        'cnvLength': args.cnvLength,
        'cnvQual': args.cnvQual,
        'cnvBinSupportRatio': args.cnvBinSupportRatio,
        'cnvCopyRatio': args.cnvCopyRatio,
        'Chromosome': args.Chromosome,
        'Classification': args.Classification,
        'is_rare': args.is_rare,
        'is_long': args.is_long,
        'is_brain': args.is_brain
    }
    
    main(args.csv_file, filters, args.output_file, args.verbose)
