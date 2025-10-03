import sys
import os
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from matplotlib import gridspec
import seaborn as sns
import gseapy as gp
import decoupler as dc
import scanpy as sc

def perform_deg_analysis(adata_subset, output_dir, layer='log1p_norm_cb', condition_column='condition', reference='WT', comparison='KO', 
                        sex_label=''):
    """
    Perform differential gene expression analysis between conditions
    
    Parameters:
    -----------
    adata_subset : AnnData object
        Subsetted data (e.g., adata_male or adata_female)
    output_dir : str
        Directory to save the results
    layer : str
        Counts layer to use for the analysis
    condition_column : str
        Column name in adata.obs that contains condition information
    reference : str
        Reference condition (e.g., 'WT')
    comparison : str
        Comparison condition (e.g., 'KO')
    sex_label : str
        Label for the sex subset ('Male' or 'Female')
    """
    
    print(f"=== Performing DEG Analysis: {sex_label} {comparison} vs {reference} ===")
    
    # Make a copy to avoid modifying original data
    adata_work = adata_subset.copy()
    
    # Ensure we have the right conditions
    conditions = adata_work.obs[condition_column].unique()
    print(f"Available conditions: {conditions}")
    
    if reference not in conditions or comparison not in conditions:
        print(f"Warning: {reference} or {comparison} not found in {condition_column}")
        return None, None
    
    # Filter to only include the two conditions of interest
    mask = adata_work.obs[condition_column].isin([reference, comparison])
    adata_filtered = adata_work[mask].copy()
    
    print(f"Cells in analysis - {reference}: {(adata_filtered.obs[condition_column] == reference).sum()}")
    print(f"Cells in analysis - {comparison}: {(adata_filtered.obs[condition_column] == comparison).sum()}")
    
    # Perform differential expression using Wilcoxon rank-sum test
    sc.tl.rank_genes_groups(
        adata_filtered,
        groupby=condition_column,
        groups=[comparison],  # Compare KO
        reference=reference,   # Against WT
        method='wilcoxon',
        pts=True,  # Calculate percentage of cells expressing the gene
        tie_correct=True,
        use_raw=False,
        layer=layer
    )
    
    # Extract results
    result_df = sc.get.rank_genes_groups_df(adata_filtered, group=comparison)
    
    # Add additional statistics
    result_df['log2FC'] = result_df['logfoldchanges']
    result_df['-log10(pval)'] = -np.log10(result_df['pvals'])
    result_df['significant'] = (result_df['pvals_adj'] < 0.05) & (np.abs(result_df['logfoldchanges']) > 0.5)
    
    print(f"Total genes analyzed: {len(result_df)}")
    print(f"Significantly upregulated genes (padj < 0.05, |log2FC| > 0.5): {(result_df['significant'] & (result_df['log2FC'] > 0)).sum()}")
    print(f"Significantly downregulated genes (padj < 0.05, |log2FC| > 0.5): {(result_df['significant'] & (result_df['log2FC'] < 0)).sum()}")
    
    return adata_filtered, result_df

def visualize_deg_results(result_df, sex_label, output_dir):
    """
    Create comprehensive visualizations for DEG results
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np
    import os
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'Differential Gene Expression Analysis - {sex_label} Cells (KO vs WT)', fontsize=16, fontweight='bold')
    
    # 1. Volcano plot
    ax1 = axes[0, 0]
    
    # Color points based on significance and fold change
    colors = []
    for _, row in result_df.iterrows():
        if row['pvals_adj'] < 0.05 and row['log2FC'] > 0.5:
            colors.append('red')  # Upregulated
        elif row['pvals_adj'] < 0.05 and row['log2FC'] < -0.5:
            colors.append('blue')  # Downregulated
        else:
            colors.append('gray')  # Not significant
    
    scatter = ax1.scatter(result_df['log2FC'], result_df['-log10(pval)'], 
                         c=colors, alpha=0.6, s=20)
    
    ax1.axvline(x=0.5, color='red', linestyle='--', alpha=0.5)
    ax1.axvline(x=-0.5, color='blue', linestyle='--', alpha=0.5)
    ax1.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    
    ax1.set_xlabel('Log2 Fold Change (KO vs WT)')
    ax1.set_ylabel('-Log10(p-value)')
    ax1.set_title('Volcano Plot')
    ax1.grid(True, alpha=0.3)
    
    # Add legend
    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label='Upregulated'),
                      Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8, label='Downregulated'),
                      Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=8, label='Not significant')]
    ax1.legend(handles=legend_elements, loc='upper right')
    
    # 2. MA plot - Use scores for x-axis since we only have pct_nz_group
    ax2 = axes[0, 1]
    
    # Use scores as the mean expression measure
    mean_expr = result_df['scores']
    xlabel_text = 'Gene Scores'
    
    ax2.scatter(mean_expr, result_df['log2FC'], c=colors, alpha=0.6, s=20)
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.5)
    ax2.axhline(y=0.5, color='red', linestyle='--', alpha=0.5)
    ax2.axhline(y=-0.5, color='blue', linestyle='--', alpha=0.5)
    ax2.set_xlabel(xlabel_text)
    ax2.set_ylabel('Log2 Fold Change (KO vs WT)')
    ax2.set_title('MA Plot')
    ax2.grid(True, alpha=0.3)
    
    # 3. Top upregulated genes
    ax3 = axes[1, 0]
    top_up = result_df[(result_df['significant']) & (result_df['log2FC'] > 0)].head(10)
    if len(top_up) > 0:
        ax3.barh(range(len(top_up)), top_up['log2FC'], color='red', alpha=0.7)
        ax3.set_yticks(range(len(top_up)))
        ax3.set_yticklabels(top_up['names'], fontsize=10)
        ax3.set_xlabel('Log2 Fold Change')
        ax3.set_title('Top 10 Upregulated Genes')
        ax3.grid(True, alpha=0.3)
    else:
        ax3.text(0.5, 0.5, 'No significantly\nupregulated genes', 
                ha='center', va='center', transform=ax3.transAxes, fontsize=12)
        ax3.set_title('Top 10 Upregulated Genes')
    
    # 4. Top downregulated genes
    ax4 = axes[1, 1]
    top_down = result_df[(result_df['significant']) & (result_df['log2FC'] < 0)].head(10)
    if len(top_down) > 0:
        ax4.barh(range(len(top_down)), top_down['log2FC'], color='blue', alpha=0.7)
        ax4.set_yticks(range(len(top_down)))
        ax4.set_yticklabels(top_down['names'], fontsize=10)
        ax4.set_xlabel('Log2 Fold Change')
        ax4.set_title('Top 10 Downregulated Genes')
        ax4.grid(True, alpha=0.3)
    else:
        ax4.text(0.5, 0.5, 'No significantly\ndownregulated genes', 
                ha='center', va='center', transform=ax4.transAxes, fontsize=12)
        ax4.set_title('Top 10 Downregulated Genes')
    
    plt.tight_layout()
    
    # Save the plot
    os.makedirs(output_dir, exist_ok=True)
    filename = f'{output_dir}/DEG_analysis_{sex_label.lower()}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Visualization saved at: {filename}")
    
    return fig

def save_deg_results(result_df, sex_label, output_dir):
    """
    Save DEG results to CSV files
    """    
    # Filter for significant results
    sig_results = result_df[result_df['significant']].copy()
    #sort the rows by pvals_adj
    sig_results = sig_results.sort_values('pvals_adj')
    # separate the rows by upregualtion or downregulation and save to different files
    sig_results_up = sig_results[sig_results['log2FC'] > 0]
    sig_results_down = sig_results[sig_results['log2FC'] < 0]
    if len(sig_results_up) > 0:
        sig_results_file_up = f'{output_dir}/DEG_results_{sex_label.lower()}_up.csv'
        sig_results_up.to_csv(sig_results_file_up, index=False)
        print(f"Significant DEG results saved at: {sig_results_file_up}")
    if len(sig_results_down) > 0:
        sig_results_file_down = f'{output_dir}/DEG_results_{sex_label.lower()}_down.csv'
        sig_results_down.to_csv(sig_results_file_down, index=False)
        print(f"Significant DEG results saved at: {sig_results_file_down}")
    
    return sig_results

def create_ranked_genelist(deg_df, log2fc_col='avg_log2FC', pval_col='p_val_adj', gene_col='gene_name', min_pval=1e-300):
    """
    Create a ranked gene list based on signed log2FC * -log10(adjusted p-value)
    """
    df = deg_df.copy()
    # Clip p-values
    df[pval_col] = df[pval_col].clip(lower=min_pval)
    # Calculate components and final score
    df['neg_log10_pval'] = -np.log10(df[pval_col])
    df['ranking_score'] = df[log2fc_col] * df['neg_log10_pval']
    # Sort by absolute ranking score
    ranked_df = df.sort_values('ranking_score', ascending=False)
    # Create final DataFrame with gene names as columns
    ranked_list = ranked_df.set_index(gene_col)[['ranking_score']].T
    return ranked_list

def simplify_term(term):
    # Convert to lowercase and split into words
    words = term.lower().replace('_', ' ').split()
    # Remove common words that don't add meaning
    stop_words = {'regulation', 'of', 'positive', 'negative', 'mediated', 'dependent', 
                 'activity', 'process', 'gobp', 'gomf', 'wp', 'reactome'}
    words = [w for w in words if w not in stop_words]
    return ' '.join(sorted(words))  # Sort words to match similar terms

def get_geneset_genes(df, msigdb_mice, output_file='geneset_genes.csv'):
    """
    Extract genes from genesets and save to CSV
    """
    # Create empty dictionary to store results
    geneset_genes = {}
    # For each pathway in the input dataframe
    for idx, row in df.iterrows():
        term = row['Term']
        # Get genes for this term from GSEA genesets (MSigDB)
        try:
            # Get leading edge genes if available
            genes = msigdb_mice[msigdb_mice['geneset'] == term]['genesymbol'].values
            geneset_genes[term] = {
                'Group': row['Group'],
                'Significance': row['-log10(pval)'],
                'Genes': ';'.join(genes)  # Join genes with semicolon for CSV
            }
        except KeyError:
            print(f"Warning: Could not find genes for {term}")
            continue
    
    # Convert to DataFrame
    result_df = pd.DataFrame.from_dict(geneset_genes, orient='index')
    result_df.index.name = 'Pathway'
    result_df.reset_index(inplace=True)
    
    # Save to CSV
    result_df.to_csv(output_file, index=False)
    
    return result_df

from pathlib import Path

def gmt_to_decoupler(pth: Path) -> pd.DataFrame:
    """Parse a gmt file to a decoupler pathway dataframe."""
    from itertools import chain, repeat

    pathways = {}

    with Path(pth).open("r") as f:
        for line in f:
            name, _, *genes = line.strip().split("\t")
            pathways[name] = genes

    return pd.DataFrame.from_records(
        chain.from_iterable(zip(repeat(k), v) for k, v in pathways.items()),
        columns=["geneset", "genesymbol"],
    )

# Calculate percentages for each condition
def get_cell_type_percentages(adata):
    wt_cells = adata[adata.obs['condition'] == 'WT'].obs['cell_type'].value_counts(normalize=True) * 100
    ko_cells = adata[adata.obs['condition'] == 'KO'].obs['cell_type'].value_counts(normalize=True) * 100
    df = pd.DataFrame({
        'Cell Type': wt_cells.index,
        'WT%': wt_cells.values.round(2),
        'ΔERCC1 KO%': ko_cells.values.round(2)
    })
    return df

def get_cell_type_percentages_by_sex(adata):
    # Calculate percentages for each combination of condition and sex
    wt_female = adata[(adata.obs['condition'] == 'WT') & (adata.obs['sex'] == 'F')].obs['cell_type'].value_counts(normalize=True) * 100
    wt_male = adata[(adata.obs['condition'] == 'WT') & (adata.obs['sex'] == 'M')].obs['cell_type'].value_counts(normalize=True) * 100
    ko_female = adata[(adata.obs['condition'] == 'KO') & (adata.obs['sex'] == 'F')].obs['cell_type'].value_counts(normalize=True) * 100
    ko_male = adata[(adata.obs['condition'] == 'KO') & (adata.obs['sex'] == 'M')].obs['cell_type'].value_counts(normalize=True) * 100
    df = pd.DataFrame({
        'Cell Type': wt_female.index,
        'WT F%': wt_female.values.round(2),
        'ΔERCC1 KO F%': ko_female.values.round(2),
        'WT M%': wt_male.values.round(2),
        'ΔERCC1 KO M%': ko_male.values.round(2)
    })
    return df