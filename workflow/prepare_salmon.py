import argparse
import pandas as pd
import os

def print_colorful_message(message, color):
    try:
        colors = {
            'red': '\033[91m',
            'green': '\033[92m',
            'yellow': '\033[93m',
            'blue': '\033[94m',
            'magenta': '\033[95m',
            'cyan': '\033[96m',
            'white': '\033[97m',
        }
        end_color = '\033[0m'
        message = message.encode('ascii', errors='replace').decode('ascii')  # 自动替换非ASCII字符
        if color not in colors:
            print(message)
        else:
            print(f"{colors[color]}{message}{end_color}")
    except Exception as e:
        print("[FALLBACK] " + str(message))
        
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
    non_numeric = eset.select_dtypes(exclude=['number']).columns
    numeric_means = eset.groupby(column_of_symbol).mean(numeric_only=True).reset_index()
    for col in non_numeric:
        if col != column_of_symbol:
            first_vals = eset.groupby(column_of_symbol)[col].first().reset_index()
            numeric_means = numeric_means.merge(first_vals, on=column_of_symbol)
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

    print("   ")
    print_colorful_message("#########################################################", "blue")
    print_colorful_message(" TPM Matrix Preparation from Salmon Output ", "cyan")
    print_colorful_message(" If you encounter issues, report at: ", "cyan")
    print_colorful_message(" https://github.com/YourLab/YourTool/issues ", "cyan")
    print_colorful_message("#########################################################", "blue")
    print(" Author: Your Name ")
    print(" Email: your.email@example.com ")
    print_colorful_message("#########################################################", "blue")
    print("   ")

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
        df = remove_duplicate_genes(df, column_of_symbol='Name')

        expr_cols = df.select_dtypes(include=['number']).columns
        print(f">>> Expression range: {df[expr_cols].min().min():.3f} to {df[expr_cols].max().max():.3f}")

        print(f">>> Saving to {output_matrix}")
        df.to_csv(output_matrix, index=False)
        print_colorful_message(">>> ✅ Done!", "green")

    except Exception as e:
        print_colorful_message(f"❌ Error occurred: {e}", "red")


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
