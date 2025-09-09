import pandas as pd
import numpy as np
import argparse
import os
from tqdm import tqdm
from importlib.resources import files
from iobrpy.utils.print_colorful_message import print_colorful_message

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Calculate Immunophenoscore (IPS)")
    parser.add_argument("--input", required=True, help="Path to expression matrix file (e.g., EXPR.txt)")
    parser.add_argument("--output", required=True, help="Output result file path (e.g., IPS_results.txt)")
    args = parser.parse_args()

    # Read input files
    _, in_ext = os.path.splitext(args.input)
    in_ext = in_ext.lower()
    if in_ext == '.csv':
        sep_in = ','
    elif in_ext in ['.tsv', '.txt']:
        sep_in = '\t'
    else:
        sep_in = ',' 
    gene_expression = pd.read_csv(args.input, sep=sep_in, index_col=0)
    # Only log2-transform if data appears unlogged (raw counts)
    if gene_expression.values.max() > 100:
        gene_expression = np.log2(gene_expression + 1)
    resource_pkg = 'iobrpy.resources'
    resource_path = files(resource_pkg).joinpath('IPS_genes.txt')
    IPS_genes = pd.read_csv(resource_path, sep="\t")

    # Check for missing genes
    missing_genes = IPS_genes[~IPS_genes['GENE'].isin(gene_expression.index)]
    if not missing_genes.empty:
        print("\nWarning: The following genes are missing or mismatched in the expression matrix:")
        print(missing_genes[['GENE', 'NAME']].to_string(index=False))

    # Calculate IPS scores
    IPS_genes = IPS_genes[IPS_genes['GENE'].isin(gene_expression.index)].reset_index(drop=True)

    sig_names = list(dict.fromkeys(IPS_genes['NAME']))
    groups   = { nm: IPS_genes.loc[IPS_genes['NAME']==nm, 'GENE'].tolist()
                 for nm in sig_names }
    w_means  = { nm: IPS_genes.loc[IPS_genes['NAME']==nm, 'WEIGHT'].mean()
                 for nm in sig_names }

    results = []
    for sample in tqdm(gene_expression.columns, desc="Processing Samples"):
        expr = gene_expression[sample]
        z = (expr.reindex(IPS_genes['GENE']) - expr.mean()) / expr.std(ddof=1)

        df_z = pd.DataFrame({
            'GENE': IPS_genes['GENE'],
            'NAME': IPS_genes['NAME'],
            'z':     z.values,
            'WEIGHT': IPS_genes['WEIGHT']
            })

        agg = df_z.groupby('NAME', sort=False).agg({'z':'mean', 'WEIGHT':'mean'}).reset_index()
        mig     = agg['z'].values
        weights = agg['WEIGHT'].values
        wg      = mig * weights

        # MHC, CP, EC, SC, AZ
        mhc = np.nanmean(wg[0:10])
        cp  = np.nanmean(wg[10:20])
        ec  = np.nanmean(wg[20:24])
        sc  = np.nanmean(wg[24:26])
        az  = mhc + cp + ec + sc

        ips = int(round(az * 10.0 / 3.0, 0))

        results.append({
            'Sample': sample,
            'MHC': mhc,
            'CP': cp,
            'EC': ec,
            'SC': sc,
            'AZ': az,
            'IPS': ips
        })

    # Prepare DataFrame
    result_df = pd.DataFrame(results)
    result_df = result_df[['Sample', 'MHC', 'EC', 'SC', 'CP', 'AZ', 'IPS']]
    result_df = result_df.rename(columns={'Sample':'ID'})
    new_cols = ['ID'] + [f"{c}_IPS" for c in result_df.columns.tolist()[1:]]
    result_df.columns = new_cols
    for col in ['MHC_IPS', 'CP_IPS', 'EC_IPS', 'SC_IPS', 'AZ_IPS']:
        result_df[col] = result_df[col].round(6)
    result_df['IPS_IPS'] = result_df['IPS_IPS'].astype(int)
    # Determine separator based on output extension
    _, ext = os.path.splitext(args.output)
    ext = ext.lower()
    if ext == '.csv':
        sep = ','
    elif ext in ['.tsv', '.txt']:
        sep = '\t'
    else:
        # Default to comma for unknown extensions
        sep = ','

    # Save results with the correct separator
    result_df.to_csv(args.output, sep=sep, index=False)
    
    print(f"\nResults saved to: {args.output} ")

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

if __name__ == "__main__":
    main()