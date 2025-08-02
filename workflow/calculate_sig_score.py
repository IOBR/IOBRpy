#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import zscore, spearmanr
import importlib_resources as pkg_resources

try:
    import gseapy as gp
except ImportError:
    gp = None

# Supported methods
signature_score_calculation_methods = {
    "pca": "pca",
    "zscore": "zscore",
    "ssgsea": "ssgsea",
    "integration": "integration"
}

# Debug flag
default_debug = False

def load_signatures(pkl_path):
    """Load signature collections from a pickle file."""
    return pd.read_pickle(pkl_path)

def preprocess_eset(eset, adjust_eset):
    """Log2-transform and drop Inf/zero‑sd genes like R's log2eset+check_eset."""
    if default_debug:
        print(f"DEBUG: Preprocess: shape={eset.shape}, adjust={adjust_eset}")
    if eset.shape[1] < 5000:
        eset = np.log2(eset + 1)
        # drop Inf rows
        finite = np.isfinite(eset).all(axis=1)
        eset = eset.loc[finite]
    if adjust_eset:
        # drop Inf and zero‑sd again
        finite = np.isfinite(eset).all(axis=1)
        eset = eset.loc[finite]
        eset = eset.loc[eset.std(axis=1) > 0]
    if default_debug:
        print(f"DEBUG: After preprocess: shape={eset.shape}")
    return eset

def filter_signatures(sig_dict, eset, min_genes):
    """Keep only signatures with at least min_genes present in eset."""
    out = {}
    for name, genes in sig_dict.items():
        present = [g for g in genes if g in eset.index]
        if default_debug:
            print(f"DEBUG: {name}: {len(present)} genes in eset")
        if len(present) >= min_genes:
            out[name] = present
    if default_debug:
        print(f"DEBUG: {len(out)} signatures retained (min_genes={min_genes})")
    return out

def sig_score_pca(eset, sig_dict, mini_gene_count, adjust_eset):
    """
    完全模拟 R sigScore(tmp, methods="PCA")：
    - 对 tmp.T（samples×genes）scale each gene (column)；
    - PCA 提取 PC1；
    - 用 sign(cor(PC1, colMeans(tmp))) 校正方向。
    """
    pdata = pd.DataFrame({'ID': eset.columns})
    # 1. log2 转换 + check_eset + optional adjust
    eset2 = preprocess_eset(eset, adjust_eset)

    # 2. 筛 signature
    min_size = max(mini_gene_count, 2)
    sigs = filter_signatures(sig_dict, eset2, min_size)

    # 3. 对每个 signature 计算
    for name, genes in sigs.items():
        # 子集：genes × samples
        tmp = eset2.loc[sorted(set(genes) & set(eset2.index))]

        # 转置：samples × genes
        mat = tmp.T

        # 按基因（列）标准化：center & scale
        mat = mat.sub(mat.mean(axis=0), axis=1)
        mat = mat.div(mat.std(axis=0, ddof=1).replace(0, np.nan), axis=1)

        # PCA（全 SVD，确定性）
        pca = PCA(n_components=1, svd_solver='full', random_state=0)
        pc1 = pca.fit_transform(mat).flatten()

        # 方向校正：sign(corr(pc1, colMeans(tmp)))
        mean_expr = tmp.mean(axis=0).values  # 每个样本上基因的平均值
        corr = np.corrcoef(pc1, mean_expr)[0, 1]
        direction = np.sign(corr) if not np.isnan(corr) else 1.0
        pdata[name] = pc1 * direction

    # 4. TME 差分
    if {'TMEscoreA_CIR','TMEscoreB_CIR'}.issubset(sigs):
        pdata['TMEscore_CIR'] = pdata['TMEscoreA_CIR'] - pdata['TMEscoreB_CIR']
    if {'TMEscoreA_plus','TMEscoreB_plus'}.issubset(sigs):
        pdata['TMEscore_plus'] = pdata['TMEscoreA_plus'] - pdata['TMEscoreB_plus']

    return pdata

def sig_score_zscore(eset, sig_dict, mini_gene_count, adjust_eset):
    """
    Python implementation matching R’s calculate_sig_score_zscore + sigScore(..., methods="zscore"):
    - Preprocess (log2, drop Inf/zero‐sd, optional adjust)
    - For each signature: subset genes, then simply take the mean across genes (colMeans)
    """
    # 1. 准备输出表
    pdata = pd.DataFrame({'ID': eset.columns})

    # 2. 预处理：log2 + drop Inf/zero‐sd, 可选 adjust
    eset2 = preprocess_eset(eset, adjust_eset)

    # 3. 筛 signatures（保留基因数 ≥ max(mini_gene_count, 2)）
    min_size = max(mini_gene_count, 2)
    sigs = filter_signatures(sig_dict, eset2, min_size)

    # 4. 对每个 signature，直接取每个样本的基因表达均值
    for name, genes in sigs.items():
        valid = sorted(set(genes) & set(eset2.index))
        mat = eset2.loc[valid]               # genes × samples
        # R sigScore(tmp, "zscore") == colMeans(tmp)
        pdata[name] = mat.mean(axis=0).values

    # 5. 组合 TME 差分列
    if {'TMEscoreA_CIR','TMEscoreB_CIR'}.issubset(sigs):
        pdata['TMEscore_CIR'] = pdata['TMEscoreA_CIR'] - pdata['TMEscoreB_CIR']
    if {'TMEscoreA_plus','TMEscoreB_plus'}.issubset(sigs):
        pdata['TMEscore_plus'] = pdata['TMEscoreA_plus'] - pdata['TMEscoreB_plus']

    return pdata

def sig_score_ssgsea(eset, sig_dict, mini_gene_count, adjust_eset, parallel_size):
    if gp is None:
        raise ImportError("gseapy required for ssGSEA")
    # Preprocess like R
    eset2 = preprocess_eset(eset, adjust_eset)
    # First filter with original threshold
    sigs = filter_signatures(sig_dict, eset2, mini_gene_count)
    # Then enforce min_size >= 5
    min_size = max(mini_gene_count, 5)
    sigs = filter_signatures(sig_dict, eset2, min_size)

    # Run ssGSEA with R-aligned parameters
    ss = gp.ssgsea(
        data=eset2,
        gene_sets=sigs,
        outdir=None,
        sample_norm_method='rank',    # rank-based kernel = Gaussian
        permutation_num=0,
        no_plot=True,
        threads=parallel_size,
        min_size=min_size,
        ssgsea_norm=True
    )

    # Pivot to samples × terms, keep Term as columns
    nes = ss.res2d.pivot(index='Term', columns='Name', values='NES').T.reset_index()
    nes.rename(columns={'Name': 'ID'}, inplace=True)
    # 5. 追加 TME 差分列，与 R 端一致
    if 'TMEscoreA_CIR' in nes.columns and 'TMEscoreB_CIR' in nes.columns:
        nes['TMEscore_CIR'] = nes['TMEscoreA_CIR'] - nes['TMEscoreB_CIR']
    if 'TMEscoreA_plus' in nes.columns and 'TMEscoreB_plus' in nes.columns:
        nes['TMEscore_plus'] = nes['TMEscoreA_plus'] - nes['TMEscoreB_plus']
    return nes

def sig_score_integration(eset, sig_dict, mini_gene_count, adjust_eset, parallel_size):
    # 1. 全局初筛：按 mini_gene_count 保留 signature
    filtered_sigs = {
        name: [g for g in genes if g in eset.index]
        for name, genes in sig_dict.items()
        if len([g for g in genes if g in eset.index]) >= mini_gene_count
    }

    # 2. 调用 PCA 和 zscore（使用预筛后的签名）
    p = sig_score_pca(eset, filtered_sigs, mini_gene_count, adjust_eset)
    p = p.set_index('ID').add_suffix('_PCA')

    z = sig_score_zscore(eset, filtered_sigs, mini_gene_count, adjust_eset)
    z = z.set_index('ID').add_suffix('_zscore')

    # 3. inline ssGSEA，避免二次硬编码过滤
    if gp is None:
        raise ImportError("gseapy required for ssGSEA")
    # 3.1 预处理一致
    eset2 = preprocess_eset(eset, adjust_eset)
    # 3.2 用 filtered_sigs 直接跑 ssGSEA（只按 mini_gene_count 过滤）
    ss = gp.ssgsea(
        data=eset2,
        gene_sets=filtered_sigs,
        outdir=None,
        sample_norm_method='rank',
        permutation_num=0,
        no_plot=True,
        threads=parallel_size,
        min_size=mini_gene_count,
        ssgsea_norm=True
    )
    nes = ss.res2d.pivot(index='Term', columns='Name', values='NES').T.reset_index()
    nes.rename(columns={'Name': 'ID'}, inplace=True)
    # 差分 TME
    if 'TMEscoreA_CIR' in nes.columns and 'TMEscoreB_CIR' in nes.columns:
        nes['TMEscore_CIR'] = nes['TMEscoreA_CIR'] - nes['TMEscoreB_CIR']
    if 'TMEscoreA_plus' in nes.columns and 'TMEscoreB_plus' in nes.columns:
        nes['TMEscore_plus'] = nes['TMEscoreA_plus'] - nes['TMEscoreB_plus']
    s = nes.set_index('ID').add_suffix('_ssGSEA')

    # 4. 横向合并三种方法
    return pd.concat([p, z, s], axis=1).reset_index()

def calculate_sig_score(eset, signature_name, method, mini_gene_count, adjust_eset, parallel_size):
    resource_pkg = 'iobrpy.resources'
    resource_path = pkg_resources.files(resource_pkg).joinpath('calculate_data.pkl')
    all_sigs = pd.read_pickle(resource_path)
    sig_dict = all_sigs.get(signature_name)
    if not isinstance(sig_dict, dict):
        raise KeyError(f"Signature '{signature_name}' must map to dict(name->genes)")
    m = method.lower()
    if m == 'pca': return sig_score_pca(eset, sig_dict, mini_gene_count, adjust_eset)
    if m == 'zscore': return sig_score_zscore(eset, sig_dict, mini_gene_count, adjust_eset)
    if m == 'ssgsea': return sig_score_ssgsea(eset, sig_dict, mini_gene_count, adjust_eset, parallel_size)
    if m == 'integration': return sig_score_integration(eset, sig_dict, mini_gene_count, adjust_eset, parallel_size)
    raise ValueError("Unknown method")

def main():
    p = argparse.ArgumentParser(description="Calculate signature scores (PCA, z-score, ssGSEA, integration).")
    p.add_argument(
        '-i', '--input',
        dest='input_matrix',
        required=True,
        help='Path to input expression matrix (CSV/TSV, genes×samples)'
    )
    p.add_argument(
        '--signature',
        required=True,
        help='Name of the signature within the pickle to use for scoring,including go_bp,go_cc,go_mf,signature_collection,signature_tme,signature_sc,signature_tumor,signature_metabolism,kegg,hallmark,reactome'
    )
    p.add_argument(
        '--method',
        default='pca',
        choices=list(signature_score_calculation_methods.values()),
        help='Scoring method to apply: "pca", "zscore", "ssgsea" or "integration"'
    )
    p.add_argument(
        '--mini_gene_count',
        type=int,
        default=3,
        help='Minimum number of genes required in a signature to be scored'
    )
    p.add_argument(
        '--adjust_eset',
        action='store_true',
        help='Whether to perform additional Inf/zero‐sd filtering after log2 transform'
    )
    p.add_argument(
        '--parallel_size',
        type=int,
        default=1,
        help='Number of threads to use for ssGSEA'
    )
    p.add_argument(
        '-o', '--output',
        dest='output_matrix',
        required=True,
        help='Path to save the resulting scores matrix'
    )
    p.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug output'
    )
    args = p.parse_args()

    if args.debug:
        default_debug = True
        print("DEBUG args:", vars(args))

    # load expression
    eset = pd.read_csv(args.input_matrix, sep=None, engine='python', index_col=0)
    if args.debug:
        print("DEBUG eset shape:", eset.shape)

    # calculate
    res = calculate_sig_score(eset, args.signature, args.method,
                              args.mini_gene_count, args.adjust_eset,
                              args.parallel_size)

    # final save (unchanged)
    # for ssgsea we want index (ID column); others index=False
    if args.method.lower() == 'ssgsea':
        res.to_csv(args.output_matrix, index=False)
    else:
        res.to_csv(args.output_matrix, index=False)

    print(f"Saved scores to {args.output_matrix}")

if __name__ == '__main__':
    main()