"""
Convert VCF files to CSV format and concatenate the results into one CSV file.

Usage:
======
    python vcf_to_cnv.py --vcf_path /path/to/vcf_files --output_dir /output_dir

    vcf_path: the directory containing VCF files.
    output_dir: the directory to save the intermediate and final CSV files.
"""

__authors__ = "Nadezhda Zhukova"
__contact__ = "nadezhda.zhukova@inserm.fr"
__copyright__ = "MIT"
__date__ = "2024-07-22"
__version__ = "1.0.0"

import logging
import os
import re
import pandas as pd
from tqdm import tqdm

# Configure logging to display information about the script's execution
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def determine_cnv_type(cnv_field):
    """
    Determine the CNV type code based on the copy number field.

    Parameters:
    cnv_field (str): Copy number field from the VCF file.

    Returns:
    int: CNV type code or -1 if there is an error in conversion.
    """
    try:
        # Attempt to convert the CNV field to an integer
        cn_value = int(cnv_field)
        return cn_value
    except ValueError:
        # Return -1 if the conversion fails
        return -1


def extract_vcf_fields(line):
    """
    Extract necessary fields from a VCF file line.

    Parameters:
    line (str): A single line from the VCF file.

    Returns:
    tuple: Extracted fields (chrom, pos, qual, info, sample_data).
    """
    fields = line.strip().split('\t')
    chrom = fields[0]
    pos = fields[1]
    qual = fields[5]
    info = fields[7]
    sample_data = fields[9]
    return chrom, pos, qual, info, sample_data


def extract_end_position(info, pos):
    """
    Extract the end position from the INFO field.

    Parameters:
    info (str): INFO field from the VCF file.
    pos (str): Position field from the VCF file.

    Returns:
    str: The end position.
    """
    end_match = re.search(r'END=(\d+)', info)
    end = end_match.group(1) if end_match else pos
    return end


def extract_cnv_type(sample_data):
    """
    Extract the CNV type from the sample data.

    Parameters:
    sample_data (str): Sample data field from the VCF file.

    Returns:
    int: CNV type code or -1 if the extraction fails.
    """
    sample_fields = sample_data.split(':')
    cnv_type_code = determine_cnv_type(
        sample_fields[2]) if len(sample_fields) > 2 else -1
    return cnv_type_code


def vcf_to_csv(vcf_filename, csv_filename):
    """
    Convert a VCF file to CSV format.

    Parameters:
    vcf_filename (str): Path to the input VCF file.
    csv_filename (str): Path to the output CSV file.
    """
    with open(vcf_filename, 'r', encoding='utf-8') as vcf_file:
        # Read all lines from the VCF file
        lines = vcf_file.readlines()

    data = []
    sample_id = None
    for line in lines:
        if line.startswith('##'):
            continue  # Skip meta-information lines
        if line.startswith('#'):
            # Extract the sample ID from the header line
            sample_id = line.strip().split('\t')[-1]
            continue

        # Extract necessary fields from the VCF line
        chrom, pos, qual, info, sample_data = extract_vcf_fields(line)

        # Extract the end position from the INFO field
        end = extract_end_position(info, pos)

        # Extract CNV type from sample data
        cnv_type_code = extract_cnv_type(sample_data)

        # Append the extracted data to the list
        data.append([sample_id, chrom, pos, end, qual, cnv_type_code])

    # Create a DataFrame from the extracted data
    df = pd.DataFrame(data,
                      columns=['id', 'chr', 'start', 'end', 'qual', 'type'])
    # Save the DataFrame to a CSV file
    df.to_csv(csv_filename, index=False)


def process_vcf_files(vcf_path, output_dir):
    """
    Process multiple VCF files and concatenate them into one CSV file.

    Parameters:
    vcf_path (str): Path to the directory containing VCF files.
    output_dir (str): Directory to save the intermediate and final CSV files.
    """
    # List all VCF files in the specified directory
    vcf_files = [os.path.join(vcf_path, f) for f in os.listdir(vcf_path)
                 if f.endswith('.vcf')]
    combined_df = pd.DataFrame()

    # Process each VCF file
    for vcf_file in tqdm(vcf_files):
        csv_file = os.path.join(
            output_dir, f"{os.path.basename(vcf_file).split('.')[0]}.csv")
        # Convert the VCF file to a CSV file
        vcf_to_csv(vcf_file, csv_file)

        # Read the CSV file into a DataFrame
        df_tmp = pd.read_csv(csv_file)
        # Concatenate the DataFrame with the combined DataFrame
        combined_df = pd.concat([combined_df, df_tmp], axis=0)

    # Reset the index of the combined DataFrame
    combined_df.reset_index(drop=True, inplace=True)
    # Save the combined DataFrame to a single CSV file
    combined_df.to_csv(os.path.join(output_dir, "all_cnv.csv"), index=False)
    logging.info('All VCF files processed and combined into %s',
                 os.path.join(output_dir, "all_cnv.csv"))


if __name__ == "__main__":
    import argparse

    # Argument parsing for command-line execution
    parser = argparse.ArgumentParser(description='Convert VCF files to CNV '
                                                 'CSV and combine.')
    parser.add_argument('--vcf_path', type=str, required=True,
                        help='Path to the directory containing VCF files.')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Directory to save the intermediate '
                             'and final CSVfiles.')

    args = parser.parse_args()

    # Create the output directory if it does not exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Process the VCF files
    process_vcf_files(args.vcf_path, args.output_dir)
