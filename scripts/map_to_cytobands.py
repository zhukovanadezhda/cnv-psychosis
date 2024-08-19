"""
Script to map Copy Number Variants to cytogenetic bands.

This script reads a VCF file and a cytogenetic band file, maps each CNV to its
corresponding cytogenetic bands, and outputs the updated VCF data with
cytogenetic band annotations.

Usage:
======
    python map_to_cytobands.py <annotated_vcf_file> <cytoband_file> <output_file>

Arguments:
==========
    annotated_vcf_file : str
        Path to the annotated VCF file (CSV format).
    cytoband_file : str
        Path to the cytogenetic band file (tab-separated format).
    output_file : str
        Path to save the updated annotated VCF file (CSV format).
"""

import sys
import logging
import pandas as pd
import pybedtools
from tqdm.auto import tqdm

def setup_logging() -> None:
    """Set up the logging configuration."""
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def cnv_to_cytoband(cnv: str, cytoband_bed: pybedtools.BedTool) -> str:
    """
    Map a CNV to its corresponding cytogenetic bands.

    Parameters:
    ===========
    cnv : str
        CNV string in the format 'chrom:start-end'.
    cytoband_bed : pybedtools.BedTool
        BedTool object containing cytogenetic bands.

    Returns:
    ========
    str
        A hyphen-separated string of cytogenetic bands.
    """
    chrom, positions = cnv.split(':')
    start, end = map(int, positions.split('-'))

    # Create a BedTool object for the CNV
    cnv_bed = pybedtools.BedTool([(chrom, start, end, "Cytoband")])

    # Intersect the CNV with the cytogenetic bands
    result = cnv_bed.intersect(cytoband_bed, wa=True, wb=True)

    # Extract the cytogenetic bands
    bands = {interval[7] for interval in result}

    return '-'.join(sorted(bands))


def main(vcf_file: str, cytoband_file: str, output_file: str) -> None:
    """
    Main function to read input files, map CNVs to cytogenetic bands, and save
    the updated VCF data.

    Parameters:
    ===========
    vcf_file : str
        Path to the VCF file (CSV format).
    cytoband_file : str
        Path to the cytogenetic band file (tab-separated format).
    output_file : str
        Path to save the updated annotated VCF file (CSV format).
    """
    try:
        logging.info("Loading VCF data from %s", vcf_file)
        all_cnv = pd.read_csv(vcf_file)

        logging.info("Loading cytogenetic band data from %s", cytoband_file)
        cytoband_df = pd.read_csv(cytoband_file, sep='\t', header=None,
                          names=['chrom', 'start', 'end', 'name', 'gieStain'])

        logging.info("Creating BedTool object from cytogenetic band data")
        cytoband_bed = pybedtools.BedTool.from_dataframe(
            cytoband_df[['chrom', 'start', 'end', 'name']]
        )

        logging.info("Mapping CNVs to cytogenetic bands")
        all_cnv['tmp'] = all_cnv.apply(
            lambda row: f"{row['Chromosome']}:{row['Start']}-{row['End']}",
            axis=1
        )

        # Apply the mapping with a progress bar
        tqdm.pandas(desc="Mapping CNVs to cytogenetic bands")
        all_cnv['Cytoband'] = all_cnv['tmp'].progress_apply(
            lambda cnv: cnv_to_cytoband(cnv, cytoband_bed)
        )

        all_cnv['Cytoband'] = all_cnv.apply(
            lambda row: f"{row['Chromosome'].split('chr')[-1]}{row['Cytoband']}",
            axis=1
        )

        all_cnv = all_cnv.drop('tmp', axis=1)

        logging.info("Saving updated annotated VCF data to %s", output_file)
        all_cnv.to_csv(output_file, index=False)

        logging.info("Processing completed successfully")

    except Exception as e:
        logging.error("Error during processing: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    setup_logging()
    if len(sys.argv) != 4:
        print("Usage: python cnv_to_cytoband.py "
              "<annotated_vcf_file> <cytoband_file> <output_file>")
        sys.exit(1)

    annotated_vcf_file = sys.argv[1]
    cytoband_file = sys.argv[2]
    output_file = sys.argv[3]

    main(annotated_vcf_file, cytoband_file, output_file)
