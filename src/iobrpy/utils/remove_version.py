# -*- coding: utf-8 -*-
"""
remove_version
--------------
Utility functions to detect and strip version suffixes from gene IDs,
e.g., 'ENSG000001.12' -> 'ENSG000001'.

Supported patterns
------------------
- Ensembl-like: ENS* prefixes with numeric body and optional ".<version>"
  Examples: ENSG000001.12, ENST000003.4, ENSMUSG00000012345.6
- (Optional) RefSeq-like: NM_/NR_/XM_/XR_/NC_/NG_/NT_/NW_/NZ_/CM_ with ".<version>"
  Example: NM_001005484.2 -> NM_001005484

Typical usage
-------------
    from iobrpy.utils.remove_version import strip_versions_in_dataframe
    df, n = strip_versions_in_dataframe(df, on_index=True, also_refseq=False)
    print(f"Stripped {n} versioned IDs")

    # If duplicates appear after stripping (e.g., ENSG000001.1 and .2 collapse):
    from iobrpy.utils.remove_version import deduplicate_after_stripping
    df = deduplicate_after_stripping(df, how='mean', on_index=True)
"""

from __future__ import annotations
import re
from typing import Iterable, Tuple, Optional, Literal
import pandas as pd

__all__ = [
    "strip_version_token",
    "strip_versions_iter",
    "strip_versions_in_index",
    "strip_versions_in_series",
    "strip_versions_in_dataframe",
    "deduplicate_after_stripping",
]

# ---------------------------------------------------------------------
# Regex patterns for versioned IDs
# ---------------------------------------------------------------------

# Ensembl-like: core prefix + digits + optional .version
# Matches: ENSG/ENST/ENSP/ENSR, ENSMUSG/ENSMUST/ENSMUSP, etc.
_ENS_RE = re.compile(r'^(ENS[A-Z]*\d{3,})(?:\.(\d+))$')

# RefSeq-like (optional): NM_/NR_/XM_/XR_/NC_/NG_/NT_/NW_/NZ_/CM_ + .version
_REFSEQ_RE = re.compile(r'^([A-Z]{2}_[0-9]+)(?:\.(\d+))$')


# ---------------------------------------------------------------------
# Core APIs
# ---------------------------------------------------------------------

def strip_version_token(gene_id: str, also_refseq: bool = False) -> str:
    """
    Strip version suffix from a single gene ID if it matches a known pattern.

    Parameters
    ----------
    gene_id : str
        Input identifier such as 'ENSG000001.12', 'NM_001005484.2', or 'TP53'.
    also_refseq : bool
        If True, also strip RefSeq versions like 'NM_... .2'.

    Returns
    -------
    str
        Cleaned ID (unchanged if no version pattern is detected).
    """
    if not isinstance(gene_id, str):
        try:
            gene_id = str(gene_id)
        except Exception:
            return gene_id

    m = _ENS_RE.match(gene_id)
    if m:
        return m.group(1)  # remove .<version>

    if also_refseq:
        r = _REFSEQ_RE.match(gene_id)
        if r:
            return r.group(1)

    return gene_id


def strip_versions_iter(ids: Iterable[str], also_refseq: bool = False) -> Tuple[list, int]:
    """
    Strip versions for an iterable of IDs.

    Returns
    -------
    (clean_list, n_stripped)
    """
    out, n = [], 0
    for x in ids:
        y = strip_version_token(x, also_refseq=also_refseq)
        if y != x:
            n += 1
        out.append(y)
    return out, n


def strip_versions_in_index(index: pd.Index, also_refseq: bool = False) -> Tuple[pd.Index, int]:
    """
    Strip versions in a pandas Index.

    Returns
    -------
    (new_index, n_stripped)
    """
    cleaned, n = strip_versions_iter(index.astype(str), also_refseq=also_refseq)
    return pd.Index(cleaned, name=index.name), n


def strip_versions_in_series(series: pd.Series, also_refseq: bool = False) -> Tuple[pd.Series, int]:
    """
    Strip versions in a pandas Series of IDs.

    Returns
    -------
    (new_series, n_stripped)
    """
    cleaned, n = strip_versions_iter(series.astype(str), also_refseq=also_refseq)
    return pd.Series(cleaned, index=series.index, name=series.name), n


def strip_versions_in_dataframe(
    df: pd.DataFrame,
    column: Optional[str] = None,
    on_index: bool = True,
    also_refseq: bool = False,
    inplace: bool = False,
) -> Tuple[pd.DataFrame, int]:
    """
    Strip versions in a DataFrame, operating on the index or a specific column.

    Parameters
    ----------
    df : DataFrame
    column : str, optional
        Column name to clean when on_index=False.
    on_index : bool
        True: operate on df.index; False: operate on 'column'.
    also_refseq : bool
        If True, also strip RefSeq versions.
    inplace : bool
        If True, modify df in place. Otherwise return a copy.

    Returns
    -------
    (clean_df, n_stripped)
    """
    if not inplace:
        df = df.copy()

    n = 0
    if on_index:
        new_index, n = strip_versions_in_index(df.index, also_refseq=also_refseq)
        df.index = new_index
    else:
        if column is None:
            raise ValueError("When on_index=False, you must provide 'column'.")
        new_col, n = strip_versions_in_series(df[column], also_refseq=also_refseq)
        df[column] = new_col

    return df, n


def deduplicate_after_stripping(
    df: pd.DataFrame,
    how: Literal['first', 'sum', 'mean'] = 'first',
    on_index: bool = True,
    column: Optional[str] = None,
) -> pd.DataFrame:
    """
    Collapse duplicates that may arise after stripping versions
    (e.g., ENSG000001.1 and ENSG000001.2 -> ENSG000001).

    Parameters
    ----------
    how : 'first' | 'sum' | 'mean'
        Aggregation strategy for duplicates; 'sum'/'mean' apply to numeric columns.
    on_index : bool
        True: deduplicate by index; False: by 'column'.
    column : str, optional
        Required when on_index=False.

    Returns
    -------
    DataFrame
    """
    if on_index:
        if not df.index.has_duplicates:
            return df
        if how == 'first':
            return df[~df.index.duplicated(keep='first')]
        return df.groupby(df.index).agg(how)
    else:
        if column is None:
            raise ValueError("When on_index=False, you must provide 'column'.")
        if not df[column].duplicated().any():
            return df
        if how == 'first':
            return df.drop_duplicates(subset=[column], keep='first').set_index(column)
        num = df.select_dtypes(include='number').copy()
        num[column] = df[column].values
        return num.groupby(column).agg(how)


# ---------------------------------------------------------------------
# Usage examples (comments only; not executed)
# ---------------------------------------------------------------------
# from iobrpy.utils.remove_version import strip_versions_in_dataframe, deduplicate_after_stripping
#
# # 1) Expression matrix with gene IDs in the index:
# df, n = strip_versions_in_dataframe(df, on_index=True, also_refseq=False)
# print(f"[remove_version] stripped {n} versioned IDs")
# df = deduplicate_after_stripping(df, how='mean', on_index=True)
#
# # 2) Table with a 'gene_id' column:
# df, n = strip_versions_in_dataframe(df, on_index=False, column='gene_id', also_refseq=True)
# df = deduplicate_after_stripping(df, how='first', on_index=False, column='gene_id')