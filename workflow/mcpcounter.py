import pandas as pd
import numpy as np
import warnings
from pathlib import Path
import argparse
from tqdm import tqdm
import importlib_resources as pkg_resources

warnings.filterwarnings("ignore")

def append_signatures(expression, markers):
    """Compute mean expression for each marker gene set."""
    result = pd.DataFrame()
    for population, genes in tqdm(markers.items(),
                                  desc="Estimating cell populations",
                                  total=len(markers)):
        common = expression.index.intersection(genes)
        if not common.empty:
            result[population] = expression.loc[common].mean(axis=0)
        else:
            print(f"⚠️ Warning: No matching genes for population {population}")
    return result

def preprocess_input(file_path):
    """Load and clean input expression matrix."""
    ext = Path(file_path).suffix.lower()
    if ext == '.csv':
        sep = ','
    elif ext in ('.tsv', '.txt'):
        sep = '\t'
    else:
        sep = None
    df = pd.read_csv(file_path, sep=sep, index_col=0, engine='python' if sep is None else None)
    print(f"\nInput data shape (before): {df.shape}")

    # Clean gene names: remove version suffixes and standardize
    df.index = df.index.str.replace(r'\..*', '', regex=True)
    df.index = df.index.str.upper().str.strip()

    # Handle duplicate genes: keep max value
    if df.index.duplicated().any():
        dup_count = df.index.duplicated(keep=False).sum()
        print(f"⚠️ Warning: {dup_count} duplicate genes found; keeping max expression per gene.")
        df = df.groupby(df.index).max()

    print(f"Input data shape (after): {df.shape}")
    print(f"First 10 gene names: {df.index[:10].tolist()}\n")
    return df

def load_signatures_from_pickle(pkl_path):
    """Load probesets and gene signatures from a local pickle file."""
    data = pd.read_pickle(pkl_path)
    # Expecting keys 'probesets' and 'genes'
    probesets = data.get('probesets')
    genes = data.get('genes')
    if probesets is None or genes is None:
        raise FileNotFoundError(f"Pickle file must contain 'probesets' and 'genes' DataFrames: {pkl_path}")
    return probesets, genes

def MCPcounter_estimate(expression, features_type):
    """Estimate cell population abundance from expression data."""
    resource_pkg = 'iobrpy.resources'
    resource_path = pkg_resources.files(resource_pkg).joinpath('mcp_data.pkl')
    probesets, genes_df = load_signatures_from_pickle(resource_path)

    if features_type == 'affy133P2_probesets':
        sig_df = probesets.copy()
        id_col = sig_df.columns[0]  # assume first column is ProbeID
        pop_col = sig_df.columns[1]

    else:
        sig_df = genes_df.copy()
        mapping = {
            'HUGO_symbols': ('HUGO symbols', 'Cell population'),
            'ENTREZ_ID': ('ENTREZID', 'Cell population'),
            'ENSEMBL_ID': ('ENSEMBL ID', 'Cell population')
        }
        id_col, pop_col = mapping[features_type]

    # Standardize gene IDs
    if id_col.lower() != 'probeid':
        sig_df[id_col] = sig_df[id_col].str.upper().str.strip()

    # Filter signatures by expression genes
    valid = sig_df[sig_df[id_col].isin(expression.index)]
    if valid.empty:
        raise ValueError("No common genes found. Check feature type and identifiers.")

    # Build marker dictionary
    markers = {}
    for pop in valid[pop_col].unique():
        genes = valid.loc[valid[pop_col] == pop, id_col].tolist()
        markers[pop] = genes
        print(f"Population '{pop}': {len(genes)} marker genes matched.")

    # Compute scores
    scores = append_signatures(expression, markers).T
    return scores

def main():
    parser = argparse.ArgumentParser(
        description='MCPcounter: cell population abundance estimation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True,
                        type=lambda x: Path(x).resolve(),
                        help='Input TSV file (rows: genes, columns: samples)')
    parser.add_argument('-f', '--features', required=True,
                        choices=['affy133P2_probesets', 'HUGO_symbols', 'ENTREZ_ID', 'ENSEMBL_ID'],
                        help='Type of gene identifiers')
    parser.add_argument('-o', '--output', default=Path('mcp_results.txt').resolve(),
                        type=lambda x: Path(x).resolve(),
                        help='Output file path')

    args = parser.parse_args()

    print("\nRunning MCPcounter with parameters:")
    print(f"Input file: {args.input}")
    print(f"Features type: {args.features}")
    print(f"Output file: {args.output}\n")

    expr = preprocess_input(args.input)
    results = MCPcounter_estimate(expr, args.features)
    results.T.to_csv(args.output, sep='\t', float_format='%.6f')

    print(f"\nAnalysis completed. Results saved to {args.output}")

if __name__ == '__main__':
    main()