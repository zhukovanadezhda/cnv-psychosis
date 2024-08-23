"""
create_comparison_graph.py

This script generates comparison boxplots and annotates significant differences 
between groups for a list of descriptors. It is designed to be flexible for 
different datasets.

Usage:
    create_comparison_graph.py --data FILE --descriptors DESCRIPTORS \\
    --groups GROUPS --output FILE

Arguments:
    --data FILE            Path to the dataset file (CSV).
    --descriptors LIST     List of descriptors/parameters to compare.
    --groups DICT          Mapping of group numbers to group names.
    --output FILE          Path to the output image file.

Example:
    python create_comparison_graph.py --data data.csv --descriptors \\
    "Age,IQ,PANSS" --groups "0:UHR-NC,1:UHR-C,2:SCZ,3:BD" --output output.png
"""

import os
import sys
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from typing import List, Dict


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate comparison graphs for different descriptors."
        )
    parser.add_argument('--data', 
                        required=True, 
                        help="Path to the dataset CSV file.")
    parser.add_argument('--descriptors', 
                        required=True, 
                        help="Comma-separated list of descriptors.")
    parser.add_argument('--groups', 
                        required=True, 
                        help=("Comma-separated list of group mappings "
                        "in the form of '0:UHR-NC,1:UHR-C,...'"))
    parser.add_argument('--output', 
                        required=True, 
                        help="Output file for the generated graph.")
    return parser.parse_args()


def parse_descriptors(descriptors_str: str) -> List[str]:
    return [desc.strip() for desc in descriptors_str.split(',')]


def parse_groups(groups_str: str) -> Dict[int, str]:
    group_mapping = {}
    for group in groups_str.split(','):
        key, value = group.split(':')
        group_mapping[int(key.strip())] = value.strip()
    return group_mapping


def create_boxplots(df: pd.DataFrame, 
                    descriptors: List[str], 
                    group_mapping: Dict[int, str], 
                    output_file: str):
    sns.set(style="darkgrid")
    phenotypes = df['Phenotype_1'].unique()
    colors = sns.color_palette("Spectral", n_colors=len(phenotypes))
    
    pairs = [
        (a, b) for i, a in enumerate(group_mapping.values()) 
        for b in list(group_mapping.values())[i+1:]
        ]

    # +2 for MADRS and Sex comparisons
    num_plots = len(descriptors) + 2  
    num_rows = (num_plots // 2) + (num_plots % 2)
    num_cols = 2

    fig, ax = plt.subplots(num_rows, num_cols, figsize=(20, num_rows * 6))

    for i, desc in enumerate(descriptors):
        sns.boxplot(x="Phenotype_1", y=desc, data=df, 
                    ax=ax[i // num_cols, i % num_cols], palette=colors)
        ax[i // num_cols, i % num_cols].set_title(desc)
        ax[i // num_cols, i % num_cols].set_xlabel("")
        ax[i // num_cols, i % num_cols].tick_params(axis='x', rotation=45)

        annotator = Annotator(ax=ax[i // num_cols, i % num_cols], 
                              pairs=pairs, data=df, x="Phenotype_1", y=desc)
        annotator.configure(test="Kruskal", hide_non_significant=True,
                            comparisons_correction="Benjamini-Hochberg", 
                            text_format='simple', loc='inside', verbose=False)
        annotator.apply_test().annotate()

    # MADRS plot
    madrs_ax = ax[-1, -2]
    sns.boxplot(x="Phenotype_1", y="MADRS", data=df, ax=madrs_ax, palette=colors)
    madrs_ax.set_title("MADRS")
    madrs_ax.set_xlabel("")
    madrs_ax.tick_params(axis='x', rotation=45)

    madrs_pairs = [(group_mapping[0], group_mapping[1]), 
                   (group_mapping[0], group_mapping[2]), 
                   (group_mapping[1], group_mapping[2])]
    annotator = Annotator(ax=madrs_ax, pairs=madrs_pairs, data=df, 
                          x="Phenotype_1", y="MADRS")
    annotator.configure(test="Kruskal", hide_non_significant=True, 
                        comparisons_correction="Benjamini-Hochberg", 
                        text_format='simple', loc='inside', verbose=False)
    annotator.apply_test().annotate()

    # Sex comparison plot
    sex_ax = ax[-1, -1]
    sns.countplot(data=df, x="Phenotype_1", hue="Sex", 
                  ax=sex_ax, palette=["#CC5675", "#557BC9"])
    sex_ax.set_title("Males vs Females")
    sex_ax.legend(["Female", "Male"])

    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()


def main():
    args = parse_arguments()

    # Load data
    df = pd.read_csv(args.data)

    # Parse descriptors and group mappings
    descriptors = parse_descriptors(args.descriptors)
    group_mapping = parse_groups(args.groups)

    # Create the boxplots
    create_boxplots(df, descriptors, group_mapping, args.output)


if __name__ == "__main__":
    main()
