"""
json_output.py - JSON output utilities for IOBRpy CLI harness.

This module provides utilities for formatting CLI output as JSON
for programmatic/agent consumption.
"""

import json
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional


class JSONEncoder(json.JSONEncoder):
    """Custom JSON encoder for IOBRpy types."""

    def default(self, obj: Any) -> Any:
        """Handle non-serializable types."""
        if isinstance(obj, Path):
            return str(obj)
        elif isinstance(obj, datetime):
            return obj.isoformat()
        elif isinstance(obj, Enum):
            return obj.value
        elif hasattr(obj, '__dict__'):
            # Dataclass-like objects
            return {k: v for k, v in obj.__dict__.items() if not k.startswith('_')}
        elif hasattr(obj, 'to_dict'):
            return obj.to_dict()
        return super().default(obj)


def to_json(
    data: Any,
    pretty: bool = True,
    sort_keys: bool = False,
    use_custom_encoder: bool = True,
) -> str:
    """
    Convert data to JSON string.

    Args:
        data: Data to serialize
        pretty: Pretty-print JSON
        sort_keys: Sort dictionary keys
        use_custom_encoder: Use custom JSONEncoder

    Returns:
        JSON string
    """
    if use_custom_encoder:
        return json.dumps(
            data,
            indent=2 if pretty else None,
            sort_keys=sort_keys,
            cls=JSONEncoder,
        )
    return json.dumps(
        data,
        indent=2 if pretty else None,
        sort_keys=sort_keys,
        default=str,
    )


def from_json(json_str: str) -> Any:
    """
    Parse JSON string to Python object.

    Args:
        json_str: JSON string

    Returns:
        Parsed Python object
    """
    return json.loads(json_str)


def success_response(
    data: Any = None,
    message: str = "Success",
    **kwargs,
) -> Dict[str, Any]:
    """
    Create a success response.

    Args:
        data: Response data
        message: Success message
        **kwargs: Additional fields

    Returns:
        Response dictionary
    """
    response = {
        'success': True,
        'message': message,
    }
    if data is not None:
        response['data'] = data
    response.update(kwargs)
    return response


def error_response(
    message: str,
    error_type: Optional[str] = None,
    **kwargs,
) -> Dict[str, Any]:
    """
    Create an error response.

    Args:
        message: Error message
        error_type: Error type/category
        **kwargs: Additional fields

    Returns:
        Response dictionary
    """
    response = {
        'success': False,
        'message': message,
    }
    if error_type:
        response['error_type'] = error_type
    response.update(kwargs)
    return response


def paginated_response(
    items: List[Any],
    total: int,
    page: int = 1,
    per_page: int = 20,
    **kwargs,
) -> Dict[str, Any]:
    """
    Create a paginated response.

    Args:
        items: Items for current page
        total: Total number of items
        page: Current page number
        per_page: Items per page
        **kwargs: Additional fields

    Returns:
        Response dictionary
    """
    return {
        'success': True,
        'data': items,
        'pagination': {
            'total': total,
            'page': page,
            'per_page': per_page,
            'total_pages': (total + per_page - 1) // per_page,
        },
        **kwargs,
    }


def format_table(
    headers: List[str],
    rows: List[List[Any]],
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Format tabular data as JSON.

    Args:
        headers: Column headers
        rows: Data rows
        title: Optional table title

    Returns:
        Formatted table dictionary
    """
    return {
        'success': True,
        'type': 'table',
        'title': title,
        'headers': headers,
        'rows': rows,
        'row_count': len(rows),
        'column_count': len(headers),
    }


def format_file_info(file_path: Path) -> Dict[str, Any]:
    """
    Format file information as JSON.

    Args:
        file_path: Path to file

    Returns:
        File information dictionary
    """
    if not file_path.exists():
        return error_response(f"File not found: {file_path}")

    stat = file_path.stat()
    return {
        'path': str(file_path),
        'name': file_path.name,
        'size': stat.st_size,
        'size_human': _human_readable_size(stat.st_size),
        'modified': datetime.fromtimestamp(stat.st_mtime).isoformat(),
        'created': datetime.fromtimestamp(stat.st_ctime).isoformat(),
        'is_file': file_path.is_file(),
        'is_dir': file_path.is_dir(),
    }


def format_directory_info(dir_path: Path, recursive: bool = False) -> Dict[str, Any]:
    """
    Format directory information as JSON.

    Args:
        dir_path: Path to directory
        recursive: Include all subdirectories

    Returns:
        Directory information dictionary
    """
    if not dir_path.exists():
        return error_response(f"Directory not found: {dir_path}")

    if not dir_path.is_dir():
        return error_response(f"Not a directory: {dir_path}")

    files = []
    dirs = []

    items = dir_path.rglob('*') if recursive else dir_path.iterdir()
    for item in items:
        if item.is_file():
            files.append({
                'name': item.name,
                'path': str(item),
                'size': item.stat().st_size,
                'size_human': _human_readable_size(item.stat().st_size),
            })
        elif item.is_dir() and item != dir_path:
            dirs.append({
                'name': item.name,
                'path': str(item),
            })

    return {
        'success': True,
        'path': str(dir_path),
        'name': dir_path.name,
        'is_dir': True,
        'file_count': len(files),
        'directory_count': len(dirs),
        'files': files,
        'directories': dirs,
    }


def _human_readable_size(size: int) -> str:
    """Convert bytes to human-readable string."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} PB"


def format_progress(
    current: int,
    total: int,
    message: Optional[str] = None,
    percentage: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Format progress information as JSON.

    Args:
        current: Current progress value
        total: Total value
        message: Optional progress message
        percentage: Optional percentage override

    Returns:
        Progress dictionary
    """
    if percentage is None and total > 0:
        percentage = (current / total) * 100
    else:
        percentage = 0

    return {
        'type': 'progress',
        'current': current,
        'total': total,
        'percentage': round(percentage, 2),
        'message': message,
    }
