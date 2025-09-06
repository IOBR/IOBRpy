import argparse
import pandas as pd
import numpy as np
import pickle
import math
from tqdm import tqdm
import os
from importlib.resources import files

def infer_sep(filepath: str) -> str:
    ext = os.path.splitext(filepath)[1].lower()
    if ext == '.csv':
        return ','
    elif ext in ('.tsv', '.txt'):
        return '\t'
    else:
        return '\t'

def output_matrix(df: pd.DataFrame, out_file: str):
    sep = infer_sep(out_file)
    print(f" Writing output matrix with {df.shape[0]} rows and {df.shape[1]} columns to {out_file}")
    df_out = df.copy()
    df_out.insert(0, 'NAME', df_out.index)
    df_out.to_csv(out_file, sep=sep, index=False, float_format='%.6f')

def filter_common_genes(input_df: pd.DataFrame, common_genes: pd.DataFrame):
    print(f" common_genes shape: {common_genes.shape}")
    print(f" input_df shape: {input_df.shape}")
    # Merge on GeneSymbol and row index
    merged = pd.merge(common_genes, input_df, left_on='GeneSymbol', right_index=True, how='inner')
    print(f" After merge, merged shape: {merged.shape}")
    if merged.empty:
        print(f" WARNING: No overlapping genes found")
    merged.index = merged['GeneSymbol']
    merged = merged.drop(columns=common_genes.columns)
    print(f" Filtered merged has {merged.shape[0]} genes")
    return merged

def estimate_score(input_df: pd.DataFrame, platform: str):
    print(" Starting estimate_score")
    resource_pkg = 'iobrpy.resources'
    txt_path = files(resource_pkg).joinpath('common_genes.txt')
    common_genes = pd.read_csv(txt_path, sep='\t', header=0, dtype=str)
    pkl_path = files(resource_pkg).joinpath('estimate_data.pkl')
    with open(pkl_path, 'rb') as f:
        data = pickle.load(f)
    SI_geneset = data['SI_geneset']
    PurityDataAffy = data['PurityDataAffy']
    merged_df = filter_common_genes(input_df, common_genes)
    print(f"Merged dataset includes {merged_df.shape[0]} genes.")
    m = merged_df.values
    gene_names = merged_df.index.tolist()
    sample_names = merged_df.columns.tolist()
    Ns = len(sample_names)
    Ng = len(gene_names)
    # Rank normalization
    ranks = np.zeros_like(m)
    for j in range(Ns):
        ranks[:, j] = pd.Series(m[:, j]).rank(method='average').values
    m_ranked = 10000 * ranks / Ng
    # Gene set enrichment
    gs = SI_geneset.iloc[:, 1:].values.tolist()
    gs_names = SI_geneset.index.tolist()
    score_matrix = np.zeros((len(gs), Ns))
    bar = tqdm(total=len(gs), desc="Gene-set enrichment")
    for i, gene_set in enumerate(gs):
        overlap = set(gene_set).intersection(gene_names)
        tqdm.write(f" Gene set '{gs_names[i]}' overlap: {len(overlap)} genes")
        if not overlap:
            score_matrix[i, :] = np.nan
            bar.update(1)
            continue
        set_indices = [gene_names.index(g) for g in overlap]
        ES_vec = []
        for j in range(Ns):
            col = m_ranked[:, j]
            order = np.argsort(-col)
            correl = np.abs(col[order]) ** 0.25
            TAG = np.isin(order, set_indices).astype(int)
            no_TAG = 1 - TAG
            sum_correl = correl[TAG==1].sum()
            P0 = no_TAG / (len(order) - len(set_indices))
            F0 = np.cumsum(P0)
            Pn = TAG * correl / sum_correl
            Fn = np.cumsum(Pn)
            RES = Fn - F0
            ES_vec.append(RES.sum())
        score_matrix[i, :] = ES_vec
        bar.update(1)
    bar.close()
    score_df = pd.DataFrame(score_matrix, index=gs_names, columns=sample_names)
    est_score = score_df.sum(axis=0)
    score_df.loc['ESTIMATEScore'] = est_score
    if platform == 'affymetrix':
        def convert(x): return math.cos(0.6049872018 + 0.0001467884 * x)
        purity = est_score.apply(convert).where(lambda x: x>=0, other=np.nan)
        score_df.loc['TumorPurity'] = purity
    return score_df

def main():
    parser = argparse.ArgumentParser(description='Estimate score calculation')
    parser.add_argument('--input', '-i', required=True, help='Input matrix file (genes x samples)')
    parser.add_argument('--platform', '-p', choices=['affymetrix','agilent','illumina'], default='affymetrix',help='Specify the platform type for the input data')
    parser.add_argument('--output', '-o', required=True, help='Output scores matrix file')
    args = parser.parse_args()

    print(f" Reading input matrix from {args.input}")
    sep_in = infer_sep(args.input)
    in_df = pd.read_csv(args.input, sep=sep_in, header=0, index_col=0)

    score_df = estimate_score(in_df, args.platform)
    # Transpose for sample x metrics
    score_df = score_df.T
    output_matrix(score_df, args.output)

if __name__ == '__main__':
    main()