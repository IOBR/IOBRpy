import argparse
import pandas as pd
import os
try:
    from tqdm.auto import tqdm
except Exception:
    def tqdm(iterable, **kwargs):
        return iterable
from iobrpy.utils.print_colorful_message import print_colorful_message
        
def remove_duplicate_genes(eset: pd.DataFrame, column_of_symbol: str = 'Name') -> pd.DataFrame:
    """
    Remove duplicate genes by averaging numeric expression values
    and taking the first non-numeric annotation for duplicates.

    Args:
        eset (pd.DataFrame): Expression matrix with 'Name' column.
        column_of_symbol (str): Column name for gene symbol or ID.

    Returns:
        pd.DataFrame: Deduplicated matrix with averaged values.
    """
    # group once
    group = eset.groupby(column_of_symbol)

    # numeric columns: compute group mean per column with progress
    numeric_cols = eset.select_dtypes(include=['number']).columns.tolist()
    if numeric_cols:
        numeric_dict = {}
        for col in tqdm(numeric_cols, desc="Aggregating numeric columns", unit="col"):
            # group[col].mean() returns a Series indexed by group keys
            numeric_dict[col] = group[col].mean()
        numeric_means = pd.DataFrame(numeric_dict)
        numeric_means.index.name = column_of_symbol
        numeric_means = numeric_means.reset_index()
    else:
        # no numeric columns: create a base frame with group keys so merges later work
        numeric_means = pd.DataFrame({column_of_symbol: list(group.groups.keys())})

    # non-numeric columns: merge first non-null value per group (show progress)
    non_numeric = [c for c in eset.select_dtypes(exclude=['number']).columns if c != column_of_symbol]
    if non_numeric:
        for col in tqdm(non_numeric, desc="Merging non-numeric columns", unit="col"):
            first_vals = group[col].first().reset_index()
            numeric_means = numeric_means.merge(first_vals, on=column_of_symbol, how='left')

    return numeric_means

def prepare_salmon_tpm(eset_path: str,
                       output_matrix: str,
                       return_feature: str = 'symbol',
                       remove_version: bool = False):
    """
    Process Salmon quantification matrix and export cleaned TPM matrix.

    Args:
        eset_path (str): Path to input Salmon expression file (.tsv or .tsv.gz).
        output_matrix (str): Output path for processed expression matrix.
        return_feature (str): Gene identifier to return (symbol, ENST, ENSG).
        remove_version (bool): Whether to remove version suffix from IDs.
    """

    try:
        print(f">>> Loading file: {eset_path}")
        df = pd.read_csv(eset_path, sep='\t', compression='infer')

        if 'Name' not in df.columns:
            raise ValueError("'Name' column is missing in the input file. Ensure Salmon output is correct.")

        print(">>> Parsing annotation from 'Name' column")
        anno = df['Name'].str.split('|', expand=True).iloc[:, :8]
        anno.columns = ['ENST', 'ENSG', 'OTTHUMG', 'OTTHUMT', 'symbol2', 'symbol', 'length', 'biotype']

        print(">>> Selecting gene feature to return")
        rf = return_feature.lower()
        if rf == 'enst':
            df['Name'] = anno['ENST']
        elif rf == 'ensg':
            df['Name'] = anno['ENSG']
        elif rf == 'symbol':
            df['Name'] = anno['symbol']
        else:
            raise ValueError(f"Invalid return_feature: {return_feature}. Choose from ENST, ENSG, symbol.")

        if remove_version:
            df['Name'] = df['Name'].str.split('.').str[0]

        print(">>> Removing duplicate genes by averaging expression values")
        print(">>> Deduplicating (progress will be shown for non-numeric columns)...")
        df = remove_duplicate_genes(df, column_of_symbol='Name')

        expr_cols = df.select_dtypes(include=['number']).columns
        print(f">>> Expression range: {df[expr_cols].min().min():.3f} to {df[expr_cols].max().max():.3f}")

        print(f">>> Saving to {output_matrix}")
        df.to_csv(output_matrix, index=False)
        print_colorful_message("Done!", "green")

        print("   ")
        print_colorful_message("#########################################################", "blue")
        print_colorful_message(" IOBRpy: Immuno-Oncology Biological Research using Python ", "cyan")
        print_colorful_message(" If you encounter any issues, please report them at ", "cyan")
        print_colorful_message(" https://github.com/IOBR/IOBRpy/issues ", "cyan")
        print_colorful_message("#########################################################", "blue")
        print(" Author: Haonan Huang, Dongqiang Zeng")
        print(" Email: interlaken@smu.edu.cn ")
        print_colorful_message("#########################################################", "blue")
        print("   ")

    except Exception as e:
        print_colorful_message(f"Error occurred: {e}", "red")


def main():
    parser = argparse.ArgumentParser(description='Prepare TPM matrix from Salmon quantification output')
    parser.add_argument('--input', '-i', dest='eset_path', type=str, required=True,
                        help='Path to input Salmon file (TSV or TSV.GZ)')
    parser.add_argument('--output', '-o', dest='output_matrix', type=str, required=True,
                        help='Path to save cleaned TPM matrix')
    parser.add_argument('--return_feature', '-r', type=str, choices=['ENST', 'ENSG', 'symbol'], default='symbol',
                        help='Which gene feature to retain (ENST, ENSG, symbol)')
    parser.add_argument('--remove_version', action='store_true',
                        help='Remove version suffix from gene IDs (e.g. ENSG000001.1 → ENSG000001)')

    args = parser.parse_args()

    prepare_salmon_tpm(eset_path=args.eset_path,
                       output_matrix=args.output_matrix,
                       return_feature=args.return_feature,
                       remove_version=args.remove_version)

if __name__ == '__main__':
    main()
