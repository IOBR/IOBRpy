"""
io.py - I/O utilities for IOBRpy CLI harness.

This module provides utilities for reading and writing various file formats
commonly used in bioinformatics workflows.
"""

import json
from pathlib import Path
from typing import Optional, Union, Dict, List, Any

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


def detect_format(file_path: Path) -> Optional[str]:
    """
    Detect file format from extension.

    Args:
        file_path: Path to file

    Returns:
        Format string ('csv', 'tsv', 'txt', 'json', 'pkl', 'xlsx') or None
    """
    suffix = file_path.suffix.lower()
    if suffix == '.csv':
        return 'csv'
    elif suffix in ('.tsv', '.txt'):
        return 'tsv'
    elif suffix == '.json':
        return 'json'
    elif suffix == '.pkl':
        return 'pkl'
    elif suffix in ('.xlsx', '.xls'):
        return 'xlsx'
    elif suffix == '.gz':
        # Handle compressed files
        return detect_format(file_path.with_suffix(''))
    return None


def read_matrix(
    file_path: Path,
    index_col: int = 0,
    sep: Optional[str] = None,
    **kwargs,
) -> Optional[Union['pd.DataFrame', Dict[str, Any], List[str]]]:
    """
    Read a matrix file (CSV/TSV/JSON).

    Args:
        file_path: Path to file
        index_col: Column to use as index
        sep: Delimiter (auto-detected if None)
        **kwargs: Additional arguments for pandas read functions

    Returns:
        DataFrame, dict, or list depending on file type
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    file_format = detect_format(file_path)

    if file_format in ('csv', 'tsv'):
        if not HAS_PANDAS:
            raise ImportError("pandas is required to read CSV/TSV files")

        if sep is None:
            sep = ',' if file_format == 'csv' else '\t'

        # Handle compressed files
        compression = None
        if file_path.suffix == '.gz':
            compression = 'gzip'

        return pd.read_csv(file_path, sep=sep, index_col=index_col, compression=compression, **kwargs)

    elif file_format == 'json':
        with open(file_path, 'r') as f:
            return json.load(f)

    elif file_format == 'pkl':
        import pickle
        with open(file_path, 'rb') as f:
            return pickle.load(f)

    elif file_format == 'xlsx':
        if not HAS_PANDAS:
            raise ImportError("pandas is required to read Excel files")
        return pd.read_excel(file_path, index_col=index_col, **kwargs)

    else:
        # Try to read as text
        with open(file_path, 'r') as f:
            return f.read()


def write_matrix(
    data: Union['pd.DataFrame', Dict[str, Any], List[str], str],
    output_path: Path,
    format: Optional[str] = None,
    sep: Optional[str] = None,
    **kwargs,
) -> bool:
    """
    Write data to a file.

    Args:
        data: Data to write
        output_path: Output file path
        format: Output format (auto-detected if None)
        sep: Delimiter for CSV/TSV
        **kwargs: Additional arguments for pandas write functions

    Returns:
        True if successful
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if format is None:
        format = detect_format(output_path)

    if isinstance(data, str):
        with open(output_path, 'w') as f:
            f.write(data)
        return True

    if HAS_PANDAS and hasattr(data, 'to_csv'):
        # DataFrame-like object
        if format == 'tsv':
            data.to_csv(output_path, sep='\t', **kwargs)
        elif format == 'json':
            data.to_json(output_path, **kwargs)
        elif format == 'xlsx':
            data.to_excel(output_path, **kwargs)
        else:
            data.to_csv(output_path, sep=sep or ',', **kwargs)
        return True

    if format == 'json':
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2, default=str)
        return True

    elif format in ('csv', 'tsv'):
        if not HAS_PANDAS:
            # Write manually
            if isinstance(data, dict):
                lines = [sep.join(data.keys())]
                if all(isinstance(v, (list, tuple)) for v in data.values()):
                    # Tabular data
                    for i in range(len(list(data.values())[0])):
                        line = sep.join(str(v[i]) for v in data.values())
                        lines.append(line)
                else:
                    lines.append(sep.join(str(v) for v in data.values()))
            else:
                lines = [str(data)]
            with open(output_path, 'w') as f:
                f.write('\n'.join(lines))
            return True
        else:
            # Convert to DataFrame
            df = pd.DataFrame(data)
            df.to_csv(output_path, sep=sep or (',' if format == 'csv' else '\t'), **kwargs)
            return True

    # Default: write as text
    with open(output_path, 'w') as f:
        f.write(str(data))
    return True


def get_file_size(file_path: Path, human_readable: bool = True) -> Union[str, int]:
    """
    Get file size.

    Args:
        file_path: Path to file
        human_readable: Return human-readable string

    Returns:
        File size in bytes or human-readable string
    """
    size = file_path.stat().st_size

    if not human_readable:
        return size

    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} PB"


def count_lines(file_path: Path) -> int:
    """
    Count lines in a file.

    Args:
        file_path: Path to file

    Returns:
        Number of lines
    """
    count = 0
    with open(file_path, 'r') as f:
        for _ in f:
            count += 1
    return count


def find_files(
    directory: Path,
    pattern: str = '*',
    recursive: bool = False,
    sort_by: Optional[str] = None,
) -> List[Path]:
    """
    Find files matching a pattern.

    Args:
        directory: Directory to search
        pattern: Glob pattern
        recursive: Search recursively
        sort_by: Sort by ('name', 'size', 'modified')

    Returns:
        List of matching file paths
    """
    directory = Path(directory)

    if recursive:
        files = list(directory.rglob(pattern))
    else:
        files = list(directory.glob(pattern))

    # Filter to files only
    files = [f for f in files if f.is_file()]

    if sort_by == 'name':
        files.sort(key=lambda x: x.name)
    elif sort_by == 'size':
        files.sort(key=lambda x: x.stat().st_size)
    elif sort_by == 'modified':
        files.sort(key=lambda x: x.stat().st_mtime, reverse=True)

    return files


def ensure_directory(path: Path) -> Path:
    """
    Ensure a directory exists.

    Args:
        path: Directory path

    Returns:
        Path object
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def validate_matrix(
    data: Union['pd.DataFrame', Any],
    require_samples: int = 1,
    require_genes: int = 1,
) -> tuple[bool, str]:
    """
    Validate a matrix data structure.

    Args:
        data: Matrix data
        require_samples: Minimum number of samples (columns)
        require_genes: Minimum number of genes (rows)

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not HAS_PANDAS or not hasattr(data, 'shape'):
        return False, "Data is not a DataFrame-like object"

    if len(data.shape) != 2:
        return False, "Matrix must be 2-dimensional"

    if data.shape[0] < require_genes:
        return False, f"Matrix has {data.shape[0]} genes, minimum {require_genes} required"

    if data.shape[1] < require_samples:
        return False, f"Matrix has {data.shape[1]} samples, minimum {require_samples} required"

    return True, ""
