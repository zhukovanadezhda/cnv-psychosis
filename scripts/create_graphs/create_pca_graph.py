"""
create_pca_graph.py

This module performs PCA on a dataset and generates visualizations including 
PCA scatter plots and biplots.

Usage:
    create_pca_graph.py --data FILE --phenotype_column COLUMN \\
        --converter_column COLUMN --output FILE

Arguments:
    --data FILE                Path to the dataset file (CSV).
    --phenotype_column COLUMN  Name of the column with phenotype data.
    --converter_column COLUMN  Name of the column with converter data.
    --output FILE              Path to the output image file.

Example:
    python create_pca_graph.py --data data.csv --phenotype_column Phenotype_1 \\
        --converter_column Phenotype_2 --output output.png
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.patches import Ellipse
import numpy as np
import argparse


def parse_arguments():
    """
    Parse command-line arguments.

    Returns:
    - args (argparse.Namespace): Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Generate PCA plots for different columns in a dataset."
    )
    parser.add_argument('--data', 
                        required=True, 
                        help="Path to the dataset CSV file.")
    parser.add_argument('--phenotype_column', 
                        required=True, 
                        help="Column name for phenotype data.")
    parser.add_argument('--converter_column', 
                        required=True, 
                        help="Column name for converter data.")
    parser.add_argument('--output', 
                        required=True, 
                        help="Output file for the generated graph.")
    return parser.parse_args()


def preprocess_data(df: pd.DataFrame, columns_to_drop: list) -> tuple:
    """
    Preprocess the data by scaling features and performing PCA.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data.
    - columns_to_drop (list): List of columns to drop before PCA.

    Returns:
    - df (pd.DataFrame): DataFrame with PCA components added.
    - features (pd.DataFrame): DataFrame with the features used for PCA.
    - pca (PCA): Fitted PCA object.
    """
    features = df.drop(columns=columns_to_drop).dropna()
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features)

    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_features)

    df = df.loc[features.index].copy()
    df['PCA1'] = pca_result[:, 0]
    df['PCA2'] = pca_result[:, 1]

    return df, features, pca


def create_color_palettes(df: pd.DataFrame, 
                          phenotype_column: str, 
                          converter_column: str) -> tuple:
    """
    Create color palettes for phenotypes and converters.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data.
    - phenotype_column (str): Column name for phenotype data.
    - converter_column (str): Column name for converter data.

    Returns:
    - phenotype_palette (dict): Color palette for different phenotypes.
    - converters_palette (dict): Color palette for converters and non-converters.
    """
    phenotype_colors = sns.color_palette(
        "Spectral", 
        n_colors=len(df[phenotype_column].unique())
        )
    phenotype_palette = dict(zip(df[phenotype_column].unique(),
                                 phenotype_colors))
    converters_palette = {"converter": "#F9E29F", "non-converter": "#F18666"}

    return phenotype_palette, converters_palette


def plot_ellipse(ax, df, group_column, palette):
    """
    Plot ellipses around groups based on PCA results.

    Parameters:
    - ax (matplotlib.axes.Axes): The axis to plot on.
    - df (pd.DataFrame): DataFrame containing PCA results.
    - group_column (str): Column name with group information.
    - palette (dict): Color palette for the groups.
    """
    for group, color in palette.items():
        subset = df[df[group_column] == group]
        cov_matrix = np.cov(subset[['PCA1', 'PCA2']].T)
        mean = subset[['PCA1', 'PCA2']].mean().values

        eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
        order = eigenvalues.argsort()[::-1]
        eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]

        angle = np.arctan2(*eigenvectors[:, 0][::-1])
        angle = np.degrees(angle)

        ellipse = Ellipse(
            xy=mean,
            width=3 * np.sqrt(eigenvalues[0]),
            height=3 * np.sqrt(eigenvalues[1]),
            angle=angle,
            color=color,
            alpha=0.3
        )

        ax.add_patch(ellipse)
        

def plot_pca_scatter(ax, df, group_column, palette, title):
    """
    Plot a scatter plot of PCA results with ellipses.

    Parameters:
    - ax (matplotlib.axes.Axes): The axis to plot on.
    - df (pd.DataFrame): DataFrame containing PCA results.
    - group_column (str): Column name with group information.
    - palette (dict): Color palette for the groups.
    - title (str): Title of the plot.
    """
    sns.scatterplot(
        x='PCA1', y='PCA2',
        hue=group_column,
        palette=palette,
        data=df,
        legend="full",
        edgecolor='w',
        s=100,
        alpha=0.8,
        ax=ax
    )
    plot_ellipse(ax, df, group_column, palette)

    ax.set_title(title, fontsize=16)
    ax.set_xlabel('PCA1', fontsize=14)
    ax.set_ylabel('PCA2', fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend(title=group_column.replace("_", " ").title(), 
              title_fontsize='13', fontsize='11')


def plot_pca_biplot(ax, pca, features):
    """
    Plot a PCA biplot showing PCA components and feature vectors.

    Parameters:
    - ax (matplotlib.axes.Axes): The axis to plot on.
    - pca (PCA): Fitted PCA object.
    - features (pd.DataFrame): DataFrame with the features used for PCA.
    """
    components = pca.components_.T[:-2]
    scaling_factor = 3

    for i in range(components.shape[0]):
        ax.arrow(0, 0, components[i, 0] * scaling_factor, components[i, 1] * scaling_factor,
                  alpha=1, head_width=0.1, head_length=0.1, color="#636A9F")
        ax.text(components[i, 0] * scaling_factor * 1.33, components[i, 1] * scaling_factor * 1.18,
                 list(features.columns)[i], ha='center', va='center', fontsize=12)

    ax.set_xlim(-1.5, 2.5)
    ax.set_ylim(-1.8, 1.5)
    ax.set_title('PCA Biplot', fontsize=16)
    ax.set_xlabel('PCA1', fontsize=14)
    ax.set_ylabel('PCA2', fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.axhline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axvline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)


def generate_pca_plots(df_metadata: pd.DataFrame, 
                       phenotype_column: str, 
                       converter_column: str, 
                       output_file: str):
    """
    Generate PCA plots and save the output to a file.

    Parameters:
    - df_metadata (pd.DataFrame): DataFrame containing the metadata and features.
    - phenotype_column (str): Column name for phenotype data.
    - converter_column (str): Column name for converter data.
    - output_file (str): Path to the output image file.
    """
    sns.set_theme(style="whitegrid")
    df_metadata, features, pca = preprocess_data(
        df_metadata, 
        columns_to_drop=[phenotype_column, converter_column]
        )
    phenotype_palette, converters_palette = create_color_palettes(
        df_metadata, 
        phenotype_column, 
        converter_column
        )

    fig = plt.figure(figsize=(20, 10))
    gs = fig.add_gridspec(2, 2)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[:, 1])

    plot_pca_scatter(ax1, df_metadata, phenotype_column, phenotype_palette,
                     'PCA of Clinical Variables (Phenotypes)')
    plot_pca_scatter(ax2, df_metadata, converter_column, converters_palette,
                     'PCA of Clinical Variables (Converters vs Non-Converters)')
    plot_pca_biplot(ax3, pca, features)

    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()


def main():
    """
    Main function to execute the script from the command line.
    """
    args = parse_arguments()

    # Load data
    df = pd.read_csv(args.data)

    # Generate PCA plots
    generate_pca_plots(df, 
                       args.phenotype_column, 
                       args.converter_column, 
                       args.output)

if __name__ == "__main__":
    main()
