import argparse
import pandas as pd
import importlib_resources as pkg_res
import pickle
from pathlib import Path
import sys as _sys
from iobrpy.workflow.prepare_salmon import prepare_salmon_tpm as prepare_salmon_tpm_main
from iobrpy.workflow.count2tpm import count2tpm as count2tpm_main
from iobrpy.workflow.anno_eset import anno_eset as anno_eset_main
from iobrpy.workflow.calculate_sig_score import calculate_sig_score as calculate_sig_score_main
from iobrpy.workflow.cibersort import cibersort as cibersort_main
from iobrpy.workflow.IPS import main as IPS_main
from iobrpy.workflow.estimate import estimate_score as estimate_score_main
from iobrpy.workflow.mcpcounter import MCPcounter_estimate as MCPcounter_estimate_main
from iobrpy.workflow.mcpcounter import preprocess_input as preprocess_input_main
from iobrpy.workflow.quantiseq import main as quantiseq_main
from iobrpy.workflow.epic import main as epic_main
from iobrpy.workflow.deside import main as deside_main
from iobrpy.workflow.tme_cluster import main as tme_cluster_main

VERSION = "0.1.2"

def main():
    parser = argparse.ArgumentParser(prog='iobrpy', description="Intratumoral Microbiome Finder Tool")
    parser.add_argument('--version', action='version', version=f'iobrpy {VERSION}')

    subparsers = parser.add_subparsers(dest='command', required=True)

    # Step 1: prepare_salmon
    p1 = subparsers.add_parser('prepare_salmon', help='Prepare Salmon data matrix')
    p1.add_argument('-i', '--input', dest='eset_path', required=True,
                    help='Path to input Salmon file (TSV or TSV.GZ)')
    p1.add_argument('-o', '--output', dest='output_matrix', required=True,
                    help='Path to save cleaned TPM matrix')
    p1.add_argument('-r', '--return_feature', dest='return_feature', choices=['ENST','ENSG','symbol'], default='symbol',
                    help='Which gene feature to retain')
    p1.add_argument('--remove_version', action='store_true',
                    help='Remove version suffix from gene IDs')

    # Step 2: count2tpm
    p2 = subparsers.add_parser('count2tpm', help='Convert count matrix to TPM')
    p2.add_argument('-i', '--input', dest='count_mat', required=True,
                    help='Path to input count matrix (CSV/TSV, genes×samples)')
    p2.add_argument('--effLength_csv', type=str,
                    help='Optional CSV with id, eff_length, and gene_symbol columns')
    p2.add_argument('--idType', choices=["Ensembl","entrez","symbol","mgi"], default="Ensembl",
                    help='Gene ID type')
    p2.add_argument('--org', choices=["hsa","mmus"], default="hsa",
                    help='Organism: hsa or mmus')
    p2.add_argument('--source', choices=["local","biomart"], default="local",
                    help='Source of gene lengths')
    p2.add_argument('--id', dest='id_col', default="id",
                    help='Column name for gene ID in effLength CSV')
    p2.add_argument('--length', dest='length_col', default="eff_length",
                    help='Column name for gene length in effLength CSV')
    p2.add_argument('--gene_symbol', dest='gene_symbol_col', default="symbol",
                    help='Column name for gene symbol in effLength CSV')
    p2.add_argument('--check_data', action='store_true',
                    help='Check and remove missing values in count matrix')
    p2.add_argument('-o', '--output', dest='output_path', required=True,
                    help='Path to save TPM matrix')

    # Step 3: anno_eset
    p3 = subparsers.add_parser('anno_eset', help='Annotate expression set and remove duplicates')
    p3.add_argument('-i', '--input', dest='input_path', required=True,
                    help='Path to input expression set')
    p3.add_argument('-o', '--output', dest='output_path', required=True,
                    help='Path to save annotated expression set')
    p3.add_argument('--annotation', required=True,
                    choices=['anno_hug133plus2','anno_rnaseq','anno_illumina'],
                    help='Annotation key to use')
    p3.add_argument('--symbol', default='symbol',
                    help='Annotation symbol column')
    p3.add_argument('--probe', default='probe_id',
                    help='Annotation probe column')
    p3.add_argument('--method', default='mean', choices=['mean','sd','sum'],
                    help='Dup handling method')

    # Step 4: calculate_sig_score
    p4 = subparsers.add_parser('calculate_sig_score', help='Calculate signature scores')
    p4.add_argument('-i', '--input', dest='input_path', required=True,
                    help='Path to input expression matrix')
    p4.add_argument('-o', '--output', dest='output_path', required=True,
                    help='Path to save signature scores')
    p4.add_argument('--signature', required=True,
                    help='Signature name to use for scoring')
    p4.add_argument('--method', dest='score_method', choices=['pca','zscore','ssgsea','integration'], default='pca',
                    help='Scoring method to apply')
    p4.add_argument('--mini_gene_count', type=int, default=3,
                    help='Minimum genes per signature')
    p4.add_argument('--adjust_eset', action='store_true',
                    help='Apply additional filtering after log2 transform')
    p4.add_argument('--parallel_size', type=int, default=1,
                    help='Threads for ssGSEA')

    # Step 5: cibersort
    p5 = subparsers.add_parser('cibersort', help='Run CIBERSORT deconvolution')
    p5.add_argument('-i', '--input', dest='input_path', required=True,
                    help='Path to mixture file (CSV or TSV)')
    p5.add_argument('--perm', type=int, default=100,
                    help='Number of permutations')
    p5.add_argument('--QN', type=lambda x: x.lower() == 'true', default=True,
                    help='Quantile normalization (True/False)')
    p5.add_argument('--absolute', type=lambda x: x.lower() == 'true', default=False,
                    help='Absolute mode (True/False)')
    p5.add_argument('--abs_method', default='sig.score',
                    choices=['sig.score','no.sumto1'],
                    help='Absolute scoring method')
    p5.add_argument('--threads', type=int, dest='threads', default=1,
                    help='Number of parallel threads')
    p5.add_argument('-o', '--output', dest='output_path', required=True,
                    help='Path to save CIBERSORT results (CSV or TSV)')

    # Step 6: IPS
    p6 = subparsers.add_parser('IPS', help='Calculate Immunophenoscore (IPS)')
    p6.add_argument('-i', '--input',  dest='input_path',  required=True,
                    help='Path to expression matrix file (e.g., EXPR.txt)')
    p6.add_argument('-o', '--output', dest='output_path', required=True,
                    help='Path to save IPS results (e.g., IPS_results.txt)')

    # Step 7: estimate
    p7 = subparsers.add_parser('estimate', help='Estimate score calculation')
    p7.add_argument('-i', '--input', dest='input_path', required=True,
                    help='Path to input matrix file (genes x samples)')
    p7.add_argument('-p', '--platform', choices=['affymetrix','agilent','illumina'], default='affymetrix',
                    help='Specify the platform type for the input data')
    p7.add_argument('-o', '--output', dest='output_path', required=True,
                    help='Path to save estimate results')

    # Step 8: mcpcounter
    p8 = subparsers.add_parser('mcpcounter', help='Run MCPcounter estimation')
    p8.add_argument('-i', '--input', dest='input_path', required=True,
                    help='Path to input expression matrix (TSV, genes×samples)')
    p8.add_argument('-f', '--features', required=True,
                    choices=['affy133P2_probesets','HUGO_symbols','ENTREZ_ID','ENSEMBL_ID'],
                    help='Type of gene identifiers')
    p8.add_argument('-o', '--output', dest='output_path', required=True,
                    help='Path to save MCPcounter results (TSV)')
    
    # Step 9: quantiseq
    p9 = subparsers.add_parser('quantiseq', help='Run quanTIseq deconvolution')
    p9.add_argument('-i', '--input', dest='input', required=True,
                    help='Path to the input mixture matrix TSV/CSV file (genes x samples)')
    p9.add_argument('-o', '--output', dest='output', required=True,
                    help='Path to save the deconvolution results TSV')
    p9.add_argument('--arrays', action='store_true',
                    help='Perform quantile normalization on array data before deconvolution')
    p9.add_argument('--signame', default='TIL10',
                    help='Name of the signature set to use (default: TIL10)')
    p9.add_argument('--tumor', action='store_true',
                    help='Remove genes with high expression in tumor samples')
    p9.add_argument('--scale_mrna', dest='mRNAscale', action='store_true',
                    help='Enable mRNA scaling; use raw signature proportions otherwise')
    p9.add_argument('--method', choices=['lsei','hampel','huber','bisquare'], default='lsei',
                    help='Deconvolution method: lsei (least squares) or robust norms')
    p9.add_argument('--rmgenes', default='unassigned',
                    help="Genes to remove: 'default', 'none', or comma-separated list")

    # Step 10: epic
    p10 = subparsers.add_parser('epic', help='Run EPIC deconvolution')
    p10.add_argument('-i', '--input',  dest='input',  required=True,
                    help='Path to the bulk expression matrix (genes×samples)')
    p10.add_argument('-o', '--output', dest='output', required=True,
                    help='Path to save EPIC cell fractions (CSV/TSV)')
    p10.add_argument('--reference', choices=['TRef','BRef','both'], default='TRef',
                    help='Which reference to use for deconvolution')

    # Step 11: deside
    p11 = subparsers.add_parser('deside', help='Run DeSide deconvolution')
    p11.add_argument('-m', '--model_dir', required=True,
                     help='Path to DeSide_model directory')
    p11.add_argument('-i', '--input', required=True,
                     help='Expression file (CSV/TSV) with rows=genes and columns=samples')
    p11.add_argument('-o', '--output', required=True,
                     help='Output CSV for predicted cell fractions')
    p11.add_argument('--exp_type', choices=['TPM','log_space','linear'], default='TPM',
                     help='Input format: TPM (log2 processed), log_space (log2 TPM+1), or linear (TPM/counts)')
    p11.add_argument('--gmt', nargs='+', default=None,
                     help='Optional GMT files for pathway masking')
    p11.add_argument('--print_info', action='store_true',
                     help='Show detailed logs during prediction')
    p11.add_argument('--add_cell_type', action='store_true',
                     help='Append predicted cell type labels')
    p11.add_argument('--scaling_by_constant', action='store_true', default=True,
                     help='Enable division-by-constant scaling')
    p11.add_argument('--scaling_by_sample', action='store_true',
                     help='Enable per-sample min–max scaling')
    p11.add_argument('--one_minus_alpha', action='store_true',
                     help='Use 1−α transformation for all cell types')
    p11.add_argument('--method_adding_pathway', choices=['add_to_end','convert'], default='add_to_end',
                     help='How to integrate pathway profiles: add_to_end or convert')
    p11.add_argument('--transpose', action='store_true', default=True,
                     help='Transpose input so that rows=samples and columns=genes')
    p11.add_argument('-r', '--result_dir', default=None,
                     help='Directory to save result plots')

    # Step 12: tme_cluster
    p12 = subparsers.add_parser('tme_cluster', help='Run TME clustering')
    p12.add_argument('--input', required=True,
                     help='Path to input file (CSV/TSV/TXT)')
    p12.add_argument('--output', required=True,
                     help='Path to save clustering results (CSV/TSV/TXT)')
    p12.add_argument('--features', default=None,
                     help="Feature columns to use, e.g. '1:22' if use cibersort(excluding the sample column)")
    p12.add_argument('--pattern', default=None,
                     help='Regex to select feature columns by name')
    p12.add_argument('--id', default=None,
                     help='Column name for sample IDs (default: first column)')
    p12.add_argument('--scale', action='store_true', dest='scale', default=True,
                     help='Enable z-score scaling (default: True)')
    p12.add_argument('--no-scale', action='store_false', dest='scale',
                     help='Disable z-score scaling')
    p12.add_argument('--min_nc', type=int, default=2,
                     help='Minimum number of clusters')
    p12.add_argument('--max_nc', type=int, default=6,
                     help='Maximum number of clusters')
    p12.add_argument('--max_iter', type=int, default=10,
                     help='Maximum number of iterations for the k-means algorithm')
    p12.add_argument('--tol', type=float, default=1e-4,
                     help='Convergence tolerance for cluster center updates')
    p12.add_argument('--print_result', action='store_true', default=False,
                     help='Print intermediate KL scores and cluster counts')
    p12.add_argument('--input_sep', default=None,
                     help='Field separator for input (auto-detect if not set)')
    p12.add_argument('--output_sep', default=None,
                     help='Field separator for output (auto-detect if not set)')

    args = parser.parse_args()

    if args.command == 'prepare_salmon':
        prepare_salmon_tpm_main(
            eset_path=args.eset_path,
            output_matrix=args.output_matrix,
            return_feature=args.return_feature,
            remove_version=args.remove_version
        )

    elif args.command == 'count2tpm':
        # load count matrix
        if args.count_mat.endswith('.gz'):
            count_mat = pd.read_csv(args.count_mat, sep='\t', index_col=0, compression='gzip')
        else:
            sep = '\t' if args.count_mat.endswith(('.tsv', '.tsv.gz')) else ','
            count_mat = pd.read_csv(args.count_mat, sep=sep, index_col=0)
        eff_df = pd.read_csv(args.effLength_csv) if args.effLength_csv else None
        tpm_df = count2tpm_main(
            count_mat=count_mat,
            anno_grch38=None,
            anno_gc_vm32=None,
            idType=args.idType,
            org=args.org,
            source=args.source,
            effLength_df=eff_df,
            id_col=args.id_col,
            length_col=args.length_col,
            gene_symbol_col=args.gene_symbol_col,
            check_data=args.check_data
        )
        tpm_df.to_csv(args.output_path)
        print(f"Saved TPM matrix to {args.output_path}")

    elif args.command == 'anno_eset':
        eset_df = pd.read_csv(
            args.input_path,
            sep=None,
            engine='python',
            index_col=0
        )
        result_df = anno_eset_main(
            eset_df=eset_df,
            annotation=args.annotation,
            symbol=args.symbol,
            probe=args.probe,
            method=args.method
        )
        result_df.index.name = ''
        result_df.to_csv(args.output_path, index_label='')
        print(f"Annotated expression set saved to {args.output_path}")

    elif args.command == 'calculate_sig_score':
        ext = Path(args.input_path).suffix.lower()
        if ext == '.csv':
            eset_df = pd.read_csv(args.input_path, sep=',', index_col=0)
        elif ext == '.txt':
            eset_df = pd.read_csv(args.input_path, sep='\t', index_col=0)
        else:
            eset_df = pd.read_csv(
                args.input_path,
                sep=None,
                engine='python',
                index_col=0
            )
        scores_df = calculate_sig_score_main(
            eset_df,
            args.signature,
            args.score_method,
            args.mini_gene_count,
            args.adjust_eset,
            args.parallel_size
        )
        scores_df.to_csv(args.output_path, index=False)
        print(f"Signature scores saved to {args.output_path}")
    elif args.command == 'cibersort':
        result_df = cibersort_main(
            input_path=args.input_path,
            perm=args.perm,
            QN=args.QN,
            absolute=args.absolute,
            abs_method=args.abs_method,            
            n_jobs=args.threads
        )
        result_df.columns = [col + '_CIBERSORT' for col in result_df.columns]
        result_df.index.name = 'ID'
        delim = ',' if args.output_path.lower().endswith('.csv') else '\t'
        result_df.to_csv(args.output_path, sep=delim, index=True)
        print(f"CIBERSORT results saved to {args.output_path}")
    elif args.command == 'IPS':
        _sys_argv_orig = _sys.argv[:]
        _sys.argv = [_sys.argv[0],
                     '--input',  args.input_path,
                     '--output', args.output_path]
        IPS_main()
        _sys.argv = _sys_argv_orig
    elif args.command == 'estimate':
        # Read input matrix
        sep = '\t' if args.input_path.lower().endswith(('.tsv', '.txt')) else ','
        in_df = pd.read_csv(args.input_path, sep=sep, index_col=0)
        # Call estimate function
        score_df = estimate_score_main(in_df, args.platform)
        # Transpose and write out
        score_df = score_df.T
        score_df.columns = [col + '_estimate' for col in score_df.columns]
        out_sep = '\t' if args.output_path.lower().endswith(('.tsv', '.txt')) else ','
        score_df.to_csv(args.output_path, sep=out_sep, index_label='ID')
        print(f"Estimate scores saved to {args.output_path}")
    elif args.command == 'mcpcounter':
        expr_df = preprocess_input_main(args.input_path)
        scores_df = MCPcounter_estimate_main(expr_df, args.features)
        out_df = scores_df.T
        out_df.columns = [col.replace(' ', '_') + '_MCPcounter' for col in out_df.columns]
        out_ext = Path(args.output_path).suffix.lower()
        out_sep = ',' if out_ext == '.csv' else '\t'
        out_df.to_csv(args.output_path, sep=out_sep, index_label='ID', float_format='%.7f')
        print(f"MCPcounter results saved to {args.output_path}")
    elif args.command == 'quantiseq':
        _sys_argv_orig = _sys.argv[:]
        _sys.argv = [
            _sys.argv[0],
            '-i', args.input,
            '-o', args.output,
            *( ['--arrays'] if args.arrays else [] ),
            '--signame', args.signame,
            *( ['--tumor'] if args.tumor else [] ),
            *( ['--scale_mrna'] if args.mRNAscale else [] ),
            '--method', args.method,
            '--rmgenes', args.rmgenes,
        ]
        quantiseq_main()
        _sys_argv = _sys_argv_orig
    elif args.command == 'epic':
        _sys_argv_orig = _sys.argv[:]
        _sys.argv = [
            _sys.argv[0],
            '-i',       args.input,
            '--reference', args.reference,
            '-o',       args.output
        ]
        epic_main()
        _sys.argv = _sys_argv_orig
    elif args.command == 'deside':
        _sys_argv_orig = _sys.argv[:]
        _sys.argv = [
            _sys_argv_orig[0],
            '--model_dir', args.model_dir,
            '--input', args.input,
            '--output', args.output,
            '--exp_type', args.exp_type,
            *(['--gmt'] + args.gmt if args.gmt else []),
            *(['--print_info'] if args.print_info else []),
            *(['--add_cell_type'] if args.add_cell_type else []),
            *(['--scaling_by_constant'] if args.scaling_by_constant else []),
            *(['--scaling_by_sample'] if args.scaling_by_sample else []),
            *(['--one_minus_alpha'] if args.one_minus_alpha else []),
            '--method_adding_pathway', args.method_adding_pathway,
            *(['--transpose'] if args.transpose else []),
            *(['--result_dir', args.result_dir] if args.result_dir else []),
        ]
        deside_main()
        _sys.argv = _sys_argv_orig
    elif args.command == 'tme_cluster':
        _sys_argv_orig = _sys.argv[:]
        # build args list for tme_cluster_main
        cli = [_sys_argv_orig[0]]
        cli += ['--input', args.input]
        cli += ['--output', args.output]
        if args.features: cli += ['--features', args.features]
        if args.pattern:  cli += ['--pattern', args.pattern]
        if args.id:       cli += ['--id', args.id]
        if args.scale:    cli += ['--scale']
        else:             cli += ['--no-scale']
        cli += ['--min_nc', str(args.min_nc)]
        cli += ['--max_nc', str(args.max_nc)]
        cli += ['--max_iter', str(args.max_iter)]
        cli += ['--tol', str(args.tol)]
        if args.print_result: cli += ['--print_result']
        if args.input_sep:  cli += ['--input_sep', args.input_sep]
        if args.output_sep: cli += ['--output_sep', args.output_sep]

        _sys.argv = cli
        tme_cluster_main()
        _sys.argv = _sys_argv_orig

if __name__ == "__main__":
    main()
