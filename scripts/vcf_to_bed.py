"""

A script that converts VCF files to BED files.

Usage:
======
    python3 vcf2bed.py <input_file> <output_file>
"""

__authors__ = "Nadezhda Zhukova"
__contact__ = "nadezhda.zhukova@inserm.fr"
__copyright__ = "MIT"
__date__ = "2024-05-01"
__version__ = "1.0.0"

import sys


if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python3 vcf2bed.py <input_file> <output_file>")
        sys.exit(1)

    # Read the input file name from argument
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Open the input and output files
    with open(
        input_file, "r", encoding="utf-8"
        ) as infile, open(
            output_file, "w", encoding="utf-8"
            ) as outfile:
        for line in infile:
            if not line.startswith("#"):
                # Split the line into columns
                columns = line.strip().split("\t")
                # Extract the columns with  chr, start, end, type
                new_columns = [
                    columns[0],
                    columns[2].split(":")[-1].split("-")[-2],
                    columns[2].split(":")[-1].split("-")[-1],
                    columns[4][1:-1]
                ]
                new_line = "\t".join(new_columns)
                outfile.write(new_line + "\n")