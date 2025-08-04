import os
import argparse
import sys
import json
import pandas as pd
from deside.utility import check_dir
from deside.decon_cf import DeSide
from deside.plot import plot_predicted_result
from deside.utility.read_file import read_gene_set
from importlib_resources import files
import matplotlib.pyplot as plt

base_font = 21

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run DeSide deconvolution on bulk RNA-seq data via iobrpy CLI"
    )
    parser.add_argument('-m','--model_dir', required=True,
                        help='Path to DeSide_model directory')
    parser.add_argument('-i','--input', required=True,
                        help='Expression file (CSV/TSV) with rows=genes and columns=samples')
    parser.add_argument('-o','--output', required=True,
                        help='Output CSV for predicted cell fractions')
    parser.add_argument('--exp_type', choices=['TPM','log_space','linear'], default='TPM',
                        help='Input format: TPM (The TPM that has undergone log2 processing), log_space (log2 TPM+1), or linear (TPM/counts)')
    parser.add_argument('--gmt', nargs='+', default=None,
                        help='Optional GMT files for pathway masking')
    parser.add_argument('--print_info', action='store_true',
                        help='Show detailed logs during prediction')
    parser.add_argument('--add_cell_type', action='store_true',
                        help='Append predicted cell type labels (classification mode)')
    parser.add_argument('--scaling_by_constant', action='store_true',default = True,
                        help='Enable division-by-constant scaling (default True in CLI)')
    parser.add_argument('--scaling_by_sample', action='store_true',
                        help='Enable per-sample min–max scaling (default off in CLI)')
    parser.add_argument('--one_minus_alpha', action='store_true',
                        help='Use 1−α transformation for all cell types')
    parser.add_argument('--method_adding_pathway',
                        choices=['add_to_end','convert'],
                        default='add_to_end',
                        help='How to integrate pathway profiles: add_to_end or convert')
    parser.add_argument('--transpose', action='store_true',default = True,
                        help='Transpose input so that rows=samples and columns=genes')
    parser.add_argument('-r','--result_dir', default=None,
                        help='Directory to save result plots')
    return parser.parse_args()

def main():
    args = parse_args()

    # Validate model files
    required = ['model_DeSide.h5','key_params.txt','genes_for_gep.txt','celltypes.txt']
    for fname in required:
        if not os.path.exists(os.path.join(args.model_dir, fname)):
            sys.stderr.write(f"Error: {fname} missing in {args.model_dir}\n")
            sys.stderr.write(
            "Please download the required DeSide_model files from:\n"
            "https://figshare.com/articles/dataset/DeSide_model/25117862/1?file=44330255\n"
            )
            sys.exit(1)

    # Load hyperparameters
    with open(os.path.join(args.model_dir, 'key_params.txt')) as f:
        params = json.load(f)
    hyper_params = params.get('hyper_params', {})

    # Prepare pathway mask
    if args.gmt:
        pathway_mask = read_gene_set(args.gmt)
    elif hyper_params.get('pathway_network'):
        resource_dir = files('iobrpy').joinpath('resources')
        gmt_files = []
        for fname in [
            'c2.cp.kegg.v2023.1.Hs.symbols.gmt',
            'c2.cp.reactome.v2023.1.Hs.symbols.gmt'
        ]:
            path = resource_dir.joinpath(fname)
            if path.exists():
                gmt_files.append(str(path))
        if not gmt_files:
            sys.stderr.write(
                "Error: pathway network enabled but default GMTs missing in resources.\n")
            sys.exit(1)
        pathway_mask = read_gene_set(gmt_files)
    else:
        pathway_mask = None

    # Ensure output/plot dirs
    plot_dir = args.result_dir or os.path.dirname(args.output)
    check_dir(plot_dir)

    # Instantiate model
    model = DeSide(model_dir=args.model_dir)

    # Load expression DataFrame
    ext = os.path.splitext(args.input)[1].lower()
    if ext == '.csv':
        exp_df = pd.read_csv(args.input, sep=',', index_col=0)
    elif ext in ['.txt', '.tsv']:
        exp_df = pd.read_csv(args.input, sep='\t', index_col=0)
    else:
        sys.stderr.write("Error: unsupported file format for input.\n")
        sys.exit(1)

    # Build prediction args using official flow: let DeSide align genes internally
    predict_args = {
        'input_file': exp_df,
        'exp_type': 'TPM' if args.exp_type in ['TPM','linear'] else 'log_space',
        'transpose': args.transpose,
        'print_info': args.print_info,
        'add_cell_type': args.add_cell_type,
        'scaling_by_constant': args.scaling_by_constant,
        'scaling_by_sample': args.scaling_by_sample,
        'one_minus_alpha': args.one_minus_alpha,
        'scaling_by_sample': False,
        'scaling_by_constant': True,
        'pathway_mask': pathway_mask,
        'hyper_params': hyper_params,
        'method_adding_pathway': args.method_adding_pathway,
        'output_file_path': args.output
    }

    # Run prediction
    model.predict(**predict_args)

    # Optional plotting
    if args.result_dir:
        plt.rcParams['font.size'] = base_font * 0.5 
        plt.rcParams['axes.titlesize']    = base_font * 0.6  
        plt.rcParams['axes.labelsize']    = base_font * 0.6  
        plt.rcParams['xtick.labelsize']   = base_font * 0.5  
        plt.rcParams['ytick.labelsize']   = base_font * 0.5  
        plt.rcParams['legend.fontsize']   = base_font * 0.5
        plot_predicted_result(
            cell_frac_result_fp=args.output,
            bulk_exp_fp=args.input,
            cancer_type=None,
            model_name="DeSide",
            result_dir=args.result_dir,
            font_scale=0.1
        )
        
    print(f"DeSide deconvolution complete. Results saved to {args.output}")

if __name__ == '__main__':
    main()