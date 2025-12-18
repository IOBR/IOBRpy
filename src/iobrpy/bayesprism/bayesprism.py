# src/iobrpy/bayesprism/bayesprism.py

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Optional
from importlib.resources import files, as_file

from . import prism, extract, process_input


DATA_PACKAGE = __package__  # "iobrpy.bayesprism"


def _load_csv_from_bp_data(filename: str, **read_kwargs) -> pd.DataFrame:
    """
    Load a CSV file from the bundled BP_data directory inside iobrpy.bayesprism.
    - For matrices (sc_dat, bk_dat): we want index_col=0 by default.
    - For label files (header=None): we do NOT set index_col, so we can use iloc[:, 0].
    """
    read_kwargs.setdefault("sep", ",")

    # Only set index_col=0 if user did NOT explicitly request header=None
    if not (read_kwargs.get("header", None) is None):
        read_kwargs.setdefault("index_col", 0)

    with as_file(files(DATA_PACKAGE) / "BP_data" / filename) as path:
        return pd.read_csv(path, **read_kwargs)

def _read_user_csv(path: Path, **read_kwargs) -> pd.DataFrame:
    """Read a user-provided CSV/TSV (optionally gzipped) with automatic sep."""

    suffixes = path.suffixes
    sep = "," if ".csv" in suffixes else "\t"
    compression = "gzip" if ".gz" in suffixes else None

    return pd.read_csv(path, sep=sep, compression=compression, **read_kwargs)


def run_bayesprism(
    bulk_path: Path,
    out_dir: Path,
    n_cores: int,
    sc_dat_path: Optional[Path] = None,
    cell_state_labels_path: Optional[Path] = None,
    cell_type_labels_path: Optional[Path] = None,
) -> None:
    """
    Core pipeline: load reference from BP_data, read bulk matrix, run BayesPrism,
    and write theta / theta_cv / Z_tumor to CSV files.

    Parameters
    ----------
    bulk_path : Path
        Path to bulk expression matrix (genes x samples).
    out_dir : Path
        Output directory where result CSVs will be written.
    n_cores : int
        Number of CPU cores passed to BayesPrism (n_cores).
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Load reference data (bundled or user-provided)
    #    Make sure you have sc_dat.csv, cell_type_labels.csv,
    #    cell_state_labels.csv under BP_data if you do not pass custom files.
    if sc_dat_path is not None:
        sc_dat = _read_user_csv(sc_dat_path, header=0, index_col=0).astype(np.int32)
    else:
        sc_dat = _load_csv_from_bp_data(
            "sc_dat.csv",
            header=0,      # first row: gene names
            index_col=0,   # first column: cell IDs
        ).astype(np.int32)

    if cell_type_labels_path is not None:
        cell_type_labels = list(
            _read_user_csv(cell_type_labels_path, header=None).iloc[:, 0]
        )
    else:
        cell_type_labels = list(
            _load_csv_from_bp_data("cell_type_labels.csv", header=None).iloc[:, 0]
        )

    if cell_state_labels_path is not None:
        cell_state_labels = list(
            _read_user_csv(cell_state_labels_path, header=None).iloc[:, 0]
        )
    else:
        cell_state_labels = list(
            _load_csv_from_bp_data("cell_state_labels.csv", header=None).iloc[:, 0]
        )

    # 2) Clean up reference genes (exactly as in tutorial)
    sc_dat_filtered = process_input.cleanup_genes(
        sc_dat,
        input_type="count.matrix",
        species="hs",
        gene_group=["Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"],
        exp_cells=5,
    )
    sc_dat_filtered_pc = process_input.select_gene_type(
        sc_dat_filtered, ["protein_coding"]
    )

    # 3) Load bulk matrix: genes x samples -> transpose to samples x genes
    suffixes = bulk_path.suffixes  # e.g. ['.txt', '.gz']
    compression = "gzip" if ".gz" in suffixes else None

    if ".csv" in suffixes:
        sep = ","
    else:
        # default to tab for .tsv/.txt or others
        sep = "\t"

    bulk_df = pd.read_csv(
        bulk_path,
        sep=sep,
        index_col=0,
        compression=compression,
    )

    # bulk: genes x samples  ->  mixture: samples x genes
    mixture = bulk_df.T.astype(np.int32)

    # 4) Run BayesPrism (key fixed to 'Malignant_cells')
    key = "Malignant_cells"

    my_prism = prism.Prism.new(
        reference=sc_dat_filtered_pc,
        input_type="count.matrix",
        cell_type_labels=cell_type_labels,
        cell_state_labels=cell_state_labels,
        key=key,
        mixture=mixture,
        outlier_cut=0.01,
        outlier_fraction=0.1,
    )

    bp_res = my_prism.run(n_cores=n_cores, update_gibbs=True)

    # 5) Extract outputs (same logic as in tutorial, but cell_name=Malignant_cells)
    theta = extract.get_fraction(bp_res, which_theta="final", state_or_type="type")
    theta_cv = bp_res.posterior_theta_f.theta_cv
    Z_tumor = extract.get_exp(bp_res, state_or_type="type", cell_name=key)

    # 6) Save as CSV files in the output directory
    theta_with_suffix = theta.add_suffix("_BayesPrism")
    theta_cv_with_suffix = theta_cv.add_suffix("_BayesPrism")

    theta_with_suffix.to_csv(out_dir / "theta.csv")
    theta_cv_with_suffix.to_csv(out_dir / "theta_cv.csv")
    
    # Z_tumor is an xarray.DataArray -> convert to pandas before saving
    if hasattr(Z_tumor, "to_pandas"):
        Z_tumor_df = Z_tumor.to_pandas()
    else:
        # Fallback: wrap as DataFrame
        Z_tumor_df = pd.DataFrame(Z_tumor)

    Z_tumor_df.to_csv(out_dir / "Z_tumor.csv")

def build_parser(parser: Optional[argparse.ArgumentParser] = None) -> argparse.ArgumentParser:
    """
    Build an ArgumentParser for the `iobrpy bayesprism` command.

    This is separated so that it can be re-used or tested.
    """
    if parser is None:
        parser = argparse.ArgumentParser(
            prog="iobrpy bayesprism",
            description="Run BayesPrism deconvolution using a bundled scRNA reference.",
        )

    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        required=True,
        help=(
            "Bulk expression matrix file (genes x samples). "
            "Rows are genes, columns are samples. CSV/TSV/.gz supported."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help=(
            "Output directory for BayesPrism results "
            "(theta.csv, theta_cv.csv, Z_tumor.csv)."
        ),
    )
    parser.add_argument(
        "--threads",
        dest="threads",
        type=int,
        default=8,
        help="Number of CPU cores used by BayesPrism (n_cores).",
    )
    parser.add_argument(
        "--sc_dat",
        dest="sc_dat",
        help=(
            "Path to a custom single-cell count matrix (genes x cells). "
            "If omitted, the bundled BP_data/sc_dat.csv is used."
        ),
    )
    parser.add_argument(
        "--cell_state_labels",
        dest="cell_state_labels",
        help=(
            "Path to a custom cell_state_labels file (one label per line). "
            "If omitted, the bundled BP_data/cell_state_labels.csv is used."
        ),
    )
    parser.add_argument(
        "--cell_type_labels",
        dest="cell_type_labels",
        help=(
            "Path to a custom cell_type_labels file (one label per line). "
            "If omitted, the bundled BP_data/cell_type_labels.csv is used."
        ),
    )

    return parser


def main(argv=None) -> None:
    """
    CLI entry point for BayesPrism.

    This is called from the top-level `iobrpy` CLI, or can be called directly as
    `python -m iobrpy.bayesprism.bayesprism ...`.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    bulk_path = Path(args.input)
    out_dir = Path(args.output)
    n_cores = int(args.threads)

    sc_dat_path = Path(args.sc_dat) if args.sc_dat else None
    cell_state_labels_path = Path(args.cell_state_labels) if args.cell_state_labels else None
    cell_type_labels_path = Path(args.cell_type_labels) if args.cell_type_labels else None

    run_bayesprism(
        bulk_path,
        out_dir,
        n_cores,
        sc_dat_path=sc_dat_path,
        cell_state_labels_path=cell_state_labels_path,
        cell_type_labels_path=cell_type_labels_path,
    )

    print("BayesPrism finished.")
    print(f"Results were saved to: {out_dir.resolve()}")