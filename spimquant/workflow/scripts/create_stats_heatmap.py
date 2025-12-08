"""Create heatmap visualizations from group statistics results.

This script generates heatmaps to visualize statistical test results across
brain regions for different metrics.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def create_stats_heatmap(stats_df, label_df, metric_columns, output_path):
    """Create a heatmap showing statistics across brain regions.
    
    Parameters
    ----------
    stats_df : pd.DataFrame
        DataFrame with group statistics results
    label_df : pd.DataFrame
        DataFrame with region labels and metadata
    metric_columns : list
        List of metric columns to visualize
    output_path : str
        Path to save the output PNG file
    """
    # Merge with label information
    merged = stats_df.merge(label_df, on="index", how="left", suffixes=("", "_label"))
    
    # Determine what kind of data we have
    has_pvalues = any(col.endswith("_pval") for col in stats_df.columns)
    has_tstats = any(col.endswith("_tstat") for col in stats_df.columns)
    
    if has_pvalues or has_tstats:
        # We have statistical test results
        create_test_results_heatmap(merged, metric_columns, output_path)
    else:
        # We have summary statistics
        create_summary_heatmap(merged, metric_columns, output_path)


def create_test_results_heatmap(data, metric_columns, output_path):
    """Create heatmap for statistical test results.
    
    Shows p-values and/or t-statistics in a heatmap format.
    """
    # Find all columns with statistical results
    pval_cols = [col for col in data.columns if col.endswith("_pval")]
    tstat_cols = [col for col in data.columns if col.endswith("_tstat")]
    
    if not pval_cols and not tstat_cols:
        # No statistical results to plot
        return
    
    # Create figure with subplots for different statistics
    n_plots = (1 if pval_cols else 0) + (1 if tstat_cols else 0)
    fig, axes = plt.subplots(1, n_plots, figsize=(8 * n_plots, 12))
    
    if n_plots == 1:
        axes = [axes]
    
    plot_idx = 0
    
    # Plot p-values
    if pval_cols:
        # Prepare data matrix
        matrix_data = []
        row_labels = []
        col_labels = []
        
        for col in pval_cols:
            metric_name = col.replace("_pval", "")
            col_labels.append(metric_name)
        
        for _, row in data.iterrows():
            row_labels.append(row.get("name", f"Region {row['index']}"))
            row_data = []
            for col in pval_cols:
                val = row[col]
                # Use -log10(p-value) for better visualization
                if pd.notna(val) and val > 0:
                    row_data.append(-np.log10(val))
                else:
                    row_data.append(0)
            matrix_data.append(row_data)
        
        matrix_data = np.array(matrix_data)
        
        # Create heatmap
        ax = axes[plot_idx]
        sns.heatmap(
            matrix_data,
            ax=ax,
            yticklabels=row_labels,
            xticklabels=col_labels,
            cmap="YlOrRd",
            cbar_kws={"label": "-log10(p-value)"},
            vmin=0,
            vmax=3,  # p=0.001 is the max
        )
        
        # Add significance threshold line
        ax.axhline(y=0, color="black", linewidth=2)
        ax.set_title("Statistical Significance (-log10 p-value)")
        ax.set_xlabel("Metric")
        ax.set_ylabel("Brain Region")
        
        plot_idx += 1
    
    # Plot t-statistics
    if tstat_cols:
        # Prepare data matrix
        matrix_data = []
        row_labels = []
        col_labels = []
        
        for col in tstat_cols:
            metric_name = col.replace("_tstat", "")
            col_labels.append(metric_name)
        
        for _, row in data.iterrows():
            row_labels.append(row.get("name", f"Region {row['index']}"))
            row_data = []
            for col in tstat_cols:
                val = row[col]
                row_data.append(val if pd.notna(val) else 0)
            matrix_data.append(row_data)
        
        matrix_data = np.array(matrix_data)
        
        # Create heatmap
        ax = axes[plot_idx]
        
        # Use diverging colormap centered at 0
        vmax = np.abs(matrix_data).max()
        sns.heatmap(
            matrix_data,
            ax=ax,
            yticklabels=row_labels,
            xticklabels=col_labels,
            cmap="RdBu_r",
            cbar_kws={"label": "t-statistic"},
            center=0,
            vmin=-vmax if vmax > 0 else -1,
            vmax=vmax if vmax > 0 else 1,
        )
        
        ax.set_title("Effect Direction and Magnitude (t-statistic)")
        ax.set_xlabel("Metric")
        ax.set_ylabel("Brain Region")
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def create_summary_heatmap(data, metric_columns, output_path):
    """Create heatmap for summary statistics (means, etc.).
    
    Shows mean values across regions for different metrics.
    """
    # Find columns with mean values
    mean_cols = [col for col in data.columns if "_mean_" in col or col.endswith("_mean")]
    
    if not mean_cols:
        # Try to use the raw metric columns
        mean_cols = [col for col in metric_columns if col in data.columns]
    
    if not mean_cols:
        # Nothing to plot
        return
    
    # Prepare data matrix
    matrix_data = []
    row_labels = []
    col_labels = mean_cols
    
    for _, row in data.iterrows():
        row_labels.append(row.get("name", f"Region {row['index']}"))
        row_data = [row[col] if pd.notna(row[col]) else 0 for col in mean_cols]
        matrix_data.append(row_data)
    
    matrix_data = np.array(matrix_data)
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(10, 12))
    
    sns.heatmap(
        matrix_data,
        ax=ax,
        yticklabels=row_labels,
        xticklabels=col_labels,
        cmap="viridis",
        cbar_kws={"label": "Value"},
    )
    
    ax.set_title("Summary Statistics Across Brain Regions")
    ax.set_xlabel("Metric")
    ax.set_ylabel("Brain Region")
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    # Load statistics results
    stats_df = pd.read_csv(snakemake.input.stats_tsv, sep="\t")
    
    if len(stats_df) == 0:
        raise ValueError("No data found in group statistics TSV file")
    
    # Load label information
    label_df = pd.read_csv(snakemake.input.label_tsv, sep="\t")
    
    # Get metric columns to visualize
    metric_columns = snakemake.params.metric_columns
    
    # Create heatmap
    create_stats_heatmap(
        stats_df,
        label_df,
        metric_columns,
        snakemake.output.heatmap_png
    )


if __name__ == "__main__":
    main()
