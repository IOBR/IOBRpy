import os
import glob
import numpy as np
import pandas as pd

# tqdm progress bar (with graceful fallback if not installed)
try:
    from tqdm import tqdm  # type: ignore
except Exception:  # pragma: no cover
    class _SimpleTqdm:
        """Minimal fallback when tqdm is unavailable."""

        def __init__(self, total=None, desc=None, unit=None, dynamic_ncols=None):
            self.total = total or 0
            self.desc = desc or "progress"
            self.n = 0

        def __enter__(self):
            print(f"[{self.desc}] total={self.total}")
            return self

        def update(self, n: int = 1):
            self.n += n
            print(f"[{self.desc}] {self.n}/{self.total}", flush=True)

        def set_postfix_str(self, s: str):
            # Simple fallback: just print the latest file/sample name
            print(f"[{self.desc}] processing {s}", flush=True)

        def __exit__(self, exc_type, exc, tb):
            pass

    def tqdm(*args, **kwargs):  # type: ignore
        return _SimpleTqdm(*args, **kwargs)


# Optional: map different possible column names to unified names
# Example:
#   '#count'  -> read_count
#   'CDR3nt'  -> CDR3_dna
#   'CDR3aa'  -> CDR3_amino_acids
#   'cid'     -> consensus_id
#   'cid_full_length' -> consensus_id_complete_vdj
COLUMN_ALIASES = {
    "read_count": ["read_count", "#count", "count", "Clones"],
    "frequency": ["frequency", "frequency(proportion of read_count)"],
    "CDR3_dna": ["CDR3_dna", "CDR3nt", "CDR3.nt"],
    "CDR3_amino_acids": ["CDR3_amino_acids", "CDR3aa", "CDR3.aa"],
    "V": ["V"],
    "D": ["D"],
    "J": ["J"],
    "C": ["C"],
    "consensus_id": ["consensus_id", "cid"],
    "consensus_id_complete_vdj": [
        "consensus_id_complete_vdj",
        "cid_full_length",
    ],
}


def standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename columns of the input DataFrame to unified names
    based on COLUMN_ALIASES mapping.
    """
    col_map = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        for col in df.columns:
            if col in aliases:
                col_map[col] = canonical
                break
    return df.rename(columns=col_map)


def gini_coefficient(x: np.ndarray) -> float:
    """
    Compute Gini coefficient for a 1D array x.
    If all values are zero, return NaN.
    """
    x = np.asarray(x, dtype=float)
    if np.all(x == 0):
        return np.nan

    x = x.flatten()
    x_sorted = np.sort(x)
    n = len(x_sorted)
    cumx = np.cumsum(x_sorted)

    # Standard formula for Gini coefficient
    gini = (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n
    return gini


def process_immune_data_batch(
    path_to_result: str,
    immdata_out_csv: str,
    immune_indices_out_csv: str,
) -> pd.DataFrame:
    """
    Process all *_report.tsv files in a directory, combine them
    and compute immune diversity indices for each sample.

    Parameters
    ----------
    path_to_result : str
        Directory containing multiple *_report.tsv files.
        Each file corresponds to one sample.
    immdata_out_csv : str
        Path to save combined clone-level data (CSV).
        The path (including directory) is provided by the caller.
    immune_indices_out_csv : str
        Path to save per-sample immune indices (CSV).
        The path (including directory) is provided by the caller.

    The script does NOT hard-code any output directory,
    so it can be imported and used from another script.

    Returns
    -------
    pd.DataFrame
        DataFrame with immune indices for each sample.
    """

    # Find all files whose name ends with "_report.tsv"
    pattern = os.path.join(path_to_result, "*_report.tsv")
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(f"No *_report.tsv files found in {path_to_result!r}")

    all_samples = []

    # Progress bar over samples instead of printing columns for each file
    with tqdm(
        total=len(files),
        desc="Immune data post-processing",
        unit="sample",
        dynamic_ncols=True,
    ) as pbar:
        for file_path in files:
            # Read each TSV file
            df = pd.read_csv(file_path, sep="\t")

            # Standardize column names to unified names
            df = standardize_columns(df)

            # Extract sample name from file name by removing suffix "_report.tsv"
            sample_name = os.path.basename(file_path)
            if sample_name.endswith("_report.tsv"):
                sample_name = sample_name[: -len("_report.tsv")]

            # Add Sample column
            df["Sample"] = sample_name

            all_samples.append(df)

            # Update progress bar
            pbar.set_postfix_str(sample_name)
            pbar.update(1)

    # Combine all samples into one DataFrame (clone-level data)
    immdata_combined = pd.concat(all_samples, ignore_index=True)

    # Save combined raw data as CSV
    immdata_combined.to_csv(immdata_out_csv, index=False)

    # Check required columns for downstream calculations
    required_cols = {"Sample", "read_count", "CDR3_dna"}
    missing = required_cols - set(immdata_combined.columns)
    if missing:
        raise ValueError(f"Missing required columns in combined data: {missing}")

    results = []

    # Compute immune diversity indices per sample
    for sample_name, group in immdata_combined.groupby("Sample"):
        # Clone counts
        cnts = group["read_count"].to_numpy(dtype=float)
        # CDR3 DNA sequences
        cdr3 = group["CDR3_dna"].astype(str)
        # CDR3 length per clone
        nc = cdr3.str.len().to_numpy(dtype=float)

        Nreads = cnts.sum()

        # If total reads <= 0, fill with NaN to avoid division by zero
        if Nreads <= 0:
            Nclones = int((cnts != 0).sum())
            rec = {
                "Sample": sample_name,
                "Nreads": Nreads,
                "Nclones": Nclones,
                "Length_CDR3": np.nan,
                "Shannon_Index": np.nan,
                "Evenness": np.nan,
                "Top_clone": np.nan,
                "Second_top_clone": np.nan,
                "Rare_clone": np.nan,
                "Second_Rare_clone": np.nan,
                "Gini": np.nan,
                "Gini_Simpson": np.nan,
            }
            results.append(rec)
            continue

        # Number of non-zero clones
        Nclones = int((cnts != 0).sum())

        # CDR3 length weighted by clone counts
        length_cdr3 = float((cnts * nc).sum() / Nreads)

        # Relative clone frequencies
        p = cnts / Nreads

        # Filter strictly positive probabilities for log
        mask = p > 0
        if mask.any():
            shannon = float(-(p[mask] * np.log(p[mask])).sum())
            evenness = float(shannon / np.log(mask.sum()))
        else:
            shannon = 0.0
            evenness = np.nan

        # Gini and Gini-Simpson
        gini = gini_coefficient(p)
        gini_simpson = float(1.0 - (p ** 2).sum())

        # Top and rare clone frequencies
        p_sorted_desc = np.sort(p)[::-1]
        p_sorted_asc = np.sort(p)

        top = float(p_sorted_desc[0]) if len(p_sorted_desc) > 0 else np.nan
        second_top = float(p_sorted_desc[1]) if len(p_sorted_desc) > 1 else np.nan
        rare = float(p_sorted_asc[0]) if len(p_sorted_asc) > 0 else np.nan
        second_rare = float(p_sorted_asc[1]) if len(p_sorted_asc) > 1 else np.nan

        rec = {
            "Sample": sample_name,
            "Nreads": Nreads,
            "Nclones": Nclones,
            "Length_CDR3": length_cdr3,
            "Shannon_Index": shannon,
            "Evenness": evenness,
            "Top_clone": top,
            "Second_top_clone": second_top,
            "Rare_clone": rare,
            "Second_Rare_clone": second_rare,
            "Gini": gini,
            "Gini_Simpson": gini_simpson,
        }
        results.append(rec)

    # Combine per-sample results and save as CSV
    immune_index_results_df = pd.DataFrame(results)
    immune_index_results_df.to_csv(immune_indices_out_csv, index=False)

    return immune_index_results_df