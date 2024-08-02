"""
A script to process VCF files.

It converts them to BED format, annotates using ClassifyCNV, and merges the
results.

Usage:
======
    python annotate_vcf.py --config config.yaml [--verbose] [--keeptmp]

Config:
=======
    The script requires a configuration file (e.g., config.yaml) with the
    following structure:

    data_path: "path/to/vcf/files"
    result_path: "path/to/result/files"
    classifycnv_path: "path/to/ClassifyCNV"
    genome_build: "hg38"
"""

__authors__ = "Nadezhda Zhukova"
__contact__ = "nadezhda.zhukova@inserm.fr"
__copyright__ = "MIT"
__date__ = "2024-07-29"
__version__ = "3.0.0"

import argparse
import csv
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import List, Tuple, Optional
import yaml
import pandas as pd
from tqdm import tqdm


def setup_logging(verbose: bool) -> None:
    """Configure the logging settings."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def read_config(config_file: str) -> dict:
    """
    Read the configuration file.

    Parameters:
    config_file (str): Path to the configuration file.

    Returns:
    dict: Configuration parameters.
    """
    with open(config_file, 'r', encoding='utf-8') as file:
        return yaml.safe_load(file)


def validate_paths(config: dict) -> None:
    """
    Validate the provided paths in the configuration.

    Parameters:
    config (dict): Configuration parameters.

    Raises:
    FileNotFoundError: If any of the paths are invalid.
    """
    for key in ['classifycnv_path', 'data_path', 'result_path']:
        if not os.path.isdir(config[key]):
            raise FileNotFoundError(f"{key.upper()} is not a valid directory.")


def count_vcf_files(data_path: str) -> int:
    """
    Count the number of VCF and VCF.GZ files in the specified directory.

    Parameters:
    data_path (str): Path to the data directory.

    Returns:
    int: Number of VCF and VCF.GZ files.
    """
    return sum(
        1 for f in os.listdir(data_path) if f.endswith(('.vcf', '.vcf.gz'))
    )


def decompress_vcf_files(data_path: str) -> None:
    """
    Decompress VCF.GZ files in the specified directory.

    Parameters:
    data_path (str): Path to the data directory.
    """
    for gz_file in Path(data_path).glob('*.vcf.gz'):
        logging.debug("Decompressing %s", gz_file)
        subprocess.run(['gzip', '-dk', str(gz_file)], check=True)


def vcf_to_bed(vcf_file: str, bed_file: str) -> None:
    """
    Convert a VCF file to BED format.

    Parameters:
    vcf_file (str): Path to the input VCF file.
    bed_file (str): Path to the output BED file.
    """
    logging.debug("Converting %s to BED format at %s", vcf_file, bed_file)
    subprocess.run(['python3', 'scripts/vcf_to_bed.py', vcf_file, bed_file],
                   check=True)


def annotate_with_classifycnv(bed_file: str, classifycnv_path: str,
                              genome_build: str) -> str:
    """
    Annotate BED file using ClassifyCNV.

    Parameters:
    bed_file (str): Path to the input BED file.
    classifycnv_path (str): Path to the ClassifyCNV directory.
    genome_build (str): Genome build version.

    Returns:
    str: Path to the results directory.
    """
    logging.debug("Annotating BED file %s with ClassifyCNV", bed_file)
    result = subprocess.run(
        ['python3', os.path.join(classifycnv_path, 'ClassifyCNV.py'),
         '--infile', bed_file, '--GenomeBuild', genome_build],
        capture_output=True, text=True, check=True
    )
    for line in result.stdout.splitlines():
        if "Results saved to" in line:
            return line.split("Results saved to ")[1]
    raise RuntimeError("Error: No results found.")


def process_vcf_file(vcf_file: Path, config: dict, paths: dict) -> None:
    """
    Process a single VCF file.

    Parameters:
    vcf_file (Path): Path to the VCF file.
    config (dict): Configuration parameters.
    paths (dict): Directory paths used for processing.
    """
    logging.debug("Processing %s.", vcf_file.name)
    bed_file = paths['bed'] / f"{vcf_file.stem}.bed"
    vcf_to_bed(vcf_file, bed_file)
    result_dir = annotate_with_classifycnv(
        bed_file, config['classifycnv_path'], config['genome_build']
    )
    annotated_file = Path(result_dir) / 'Scoresheet.txt'
    shutil.copy(
        annotated_file,
        paths['annotations'] / f"annotated_{vcf_file.stem}.txt"
        )
    vcf_to_csv(vcf_file, paths['vcf'] / f"{vcf_file.stem}.csv")


def extract_vcf_fields(line: str) -> tuple:
    """
    Extract necessary fields from a VCF file line.

    Parameters:
    line (str): A single line from the VCF file.

    Returns:
    tuple: Extracted fields (chrom, pos, qual, info, sample_data).
    """
    fields = line.strip().split('\t')
    return fields[0], fields[1], fields[5], fields[6], fields[7], fields[9]


def extract_end_position(info: str, pos: str) -> str:
    """
    Extract the end position from the INFO field.

    Parameters:
    info (str): INFO field from the VCF file.
    pos (str): Position field from the VCF file.

    Returns:
    str: The end position.
    """
    match = re.search(r'END=(\d+)', info)
    return match.group(1) if match else pos


def extract_cnv_type(sample_data: str) -> int:
    """
    Extract the CNV type from the sample data.

    Parameters:
    sample_data (str): Sample data field from the VCF file.

    Returns:
    int: CNV type code or -1 if the extraction fails.
    """
    fields = sample_data.split(':')
    return int(fields[2]) if len(fields) > 2 else -1


def vcf_to_csv(vcf_filename: str, csv_filename: str) -> None:
    """
    Convert a VCF file to CSV format including quality and copy number scores.

    Parameters:
    vcf_filename (str): Path to the input VCF file.
    csv_filename (str): Path to the output CSV file.
    """
    data = []
    sample_id = None
    with open(vcf_filename, 'r', encoding='utf-8') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    sample_id = line.strip().split('\t')[-1]
                continue
            chrom, pos, qual, filter, info, sample_data = extract_vcf_fields(line)
            end = extract_end_position(info, pos)
            cnv_type_code = extract_cnv_type(sample_data)
            data.append([sample_id,
                         chrom,
                         pos,
                         end, 
                         qual,
                         filter,
                         cnv_type_code])
    df = pd.DataFrame(
        data,
        columns=['id', 'chr', 'start', 'end', 'qual', 'filter', 'type']
    )
    df.to_csv(csv_filename, index=False)


def merge_annotations(annotations_path: Path,
                      vcf_path: Path,
                      output_file: Path) -> None:
    """
    Merge all annotation files into a single CSV file.

    Parameters:
    annotations_path (Path): Directory containing the annotation files.
    vcf_path (Path): Directory containing the VCF CSV files.
    output_file (Path): Path to the output CSV file.
    """
    headers, rows = process_annotation_files(annotations_path, vcf_path)
    write_output_file(output_file, headers, rows)


def process_annotation_files(
    annotations_path: Path,
    vcf_path: Path
) -> Tuple[Optional[List[str]], List[List[str]]]:
    """
    Process annotation files and collect headers and rows for the final CSV.

    Parameters:
    annotations_path (Path): Directory containing the annotation files.
    vcf_path (Path): Directory containing the VCF CSV files.

    Returns:
    Tuple[Optional[List[str]], List[List[str]]]: Headers and rows for the CSV.
    """
    rows = []
    headers = None

    for annotation_file in Path(annotations_path).glob('annotated_*.txt'):
        file_headers, file_rows = process_single_annotation_file(
            annotation_file,
            vcf_path
        )
        if headers is None:
            headers = file_headers
        rows.extend(file_rows)

    return headers, rows


def process_single_annotation_file(
    annotation_file: Path,
    vcf_path: Path) -> Tuple[List[str], List[List[str]]]:
    """
    Process a single annotation file and collect headers and rows.

    Parameters:
    annotation_file (Path): Path to the annotation file.
    vcf_path (Path): Directory containing the VCF CSV files.

    Returns:
    Tuple[List[str], List[List[str]]]: Headers and rows for the CSV file.
    """
    rows = []
    with open(annotation_file, 'r', encoding="utf-8") as file:
        reader = csv.reader(file, delimiter='\t')
        annotations = list(reader)
        headers = ["ID"] + annotations[0][1:7] + annotations[0][-2:]

        file_id = annotation_file.stem.replace('annotated_', '')
        vcf_file_path = vcf_path / f'{file_id}.csv'
        vcf_data = (
            pd.read_csv(vcf_file_path)
            if vcf_file_path.exists()
            else pd.DataFrame()
        )

        for idx, row in enumerate(annotations[1:]):
            vcf_entry = (
                vcf_data.iloc[idx]
                if not vcf_data.empty
                else {}
            )
            qual = (
                vcf_entry.get('qual', '')
                if 'qual' in vcf_entry
                else ''
            )
            filter = (
                vcf_entry.get('filter', '')
                if 'filter' in vcf_entry
                else ''
            )
            cnv_type = (
                vcf_entry.get('type', '')
                if 'type' in vcf_entry
                else ''
            )
            rows.append(
                [file_id.split(".")[0]] + row[1:7] + row[-2:] + [qual, filter, cnv_type]
            )

    return headers, rows


def write_output_file(output_file: Path,
                      headers: List[str],
                      rows: List[List[str]]) -> None:
    """
    Write the headers and rows to the output CSV file.

    Parameters:
    output_file (Path): Path to the output CSV file.
    headers (List[str]): Headers for the CSV file.
    rows (List[List[str]]): Rows for the CSV file.
    """
    with open(output_file, 'w', newline='', encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(headers + ['Quality', 'Filter', 'Copy_number'])
        writer.writerows(rows)


def clean_up_directories(paths: dict, keep_dirs: bool) -> None:
    """
    Clean up created directories.

    Parameters:
    paths (dict): Directory paths to be cleaned up.
    keep_dirs (bool): Whether to keep directories.
    """
    if not keep_dirs:
        for path in paths.values():
            if path.is_dir():
                shutil.rmtree(path)
        logging.info("Cleanup completed.")


def delete_classifycnv_results() -> None:
    """Delete the ClassifyCNV_results directory if it exists."""
    classifycnv_results_path = Path.cwd() / 'ClassifyCNV_results'
    if classifycnv_results_path.is_dir():
        shutil.rmtree(classifycnv_results_path)
        logging.debug("ClassifyCNV_results directory deleted.")


def main(config_file: str, verbose: bool, keeptmp: bool) -> None:
    """
    Execute the VCF processing workflow.

    Parameters:
    config_file (str): Path to the configuration file.
    verbose (bool): Verbosity flag for logging.
    keeptmp (bool): Whether to keep temporary directories after processing.
    """
    setup_logging(verbose)
    config = read_config(config_file)
    validate_paths(config)
    paths = {
        'vcf': Path(config['data_path']) / 'VCF',
        'bed': Path(config['data_path']) / 'BED',
        'annotations': Path(config['data_path']) / 'ANNOTATIONS'
    }
    # Create all the temporary directories
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)

    decompress_vcf_files(config['data_path'])

    vcf_files = (
        list(paths['vcf'].glob('*.vcf')) +
        list(Path(config['data_path']).glob('*.vcf'))
    )

    # Making complete logs if verbose
    if verbose:
        vcf_files_count = count_vcf_files(config['data_path'])
        processed_files = 0
        for vcf_file in vcf_files:
            process_vcf_file(vcf_file, config, paths)
            processed_files += 1
            logging.info(
                "Processed %d of %d files.",
                processed_files,
                vcf_files_count
            )
    # Making a progress bar if not verbose
    else:
        for vcf_file in tqdm(vcf_files, desc="Processing VCF files"):
            process_vcf_file(vcf_file, config, paths)

    merge_annotations(
        paths['annotations'],
        paths['vcf'],
        Path(config['result_path']) / 'merged_annotations.csv'
    )

    if not keeptmp:
        clean_up_directories(paths, keeptmp)

    delete_classifycnv_results()
    logging.info("Processing completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='VCF Processing Script')
    parser.add_argument('--config', type=str, required=True,
                        help='Path to the configuration file.')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose mode for logging.')
    parser.add_argument('--keeptmp', action='store_true',
                        help='Keep created directories after processing')
    args = parser.parse_args()
    main(args.config, args.verbose, args.keeptmp)
