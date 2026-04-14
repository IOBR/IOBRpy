"""
project.py - Project management for IOBRpy CLI harness.

This module provides project creation, configuration, and lifecycle management
for IOBRpy TME analysis projects.
"""

import json
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, List, Any


@dataclass
class ProjectConfig:
    """Configuration for an IOBRpy project."""

    name: str
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    description: Optional[str] = None
    workspace_dir: Optional[Path] = None
    input_type: Optional[str] = None  # 'fastq', 'tpm', 'counts'
    organism: str = 'hsa'  # 'hsa' or 'mmus'
    mode: Optional[str] = None  # 'salmon' or 'star' (for fastq input)
    threads: int = 8
    batch_size: int = 1

    # Paths
    fastq_dir: Optional[Path] = None
    tpm_matrix: Optional[Path] = None
    count_matrix: Optional[Path] = None
    salmon_index: Optional[Path] = None
    star_index: Optional[Path] = None

    # Analysis options
    deconvolution_methods: List[str] = field(default_factory=lambda: [
        'cibersort', 'IPS', 'estimate', 'mcpcounter', 'quantiseq', 'epic'
    ])
    signature_groups: List[str] = field(default_factory=lambda: ['all'])
    lr_cancer_type: str = 'pancan'

    # Custom settings
    custom_settings: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary, handling Path objects."""
        d = asdict(self)
        # Convert Path objects to strings for JSON serialization
        for key, value in d.items():
            if isinstance(value, Path):
                d[key] = str(value)
            elif isinstance(value, list):
                d[key] = [str(v) if isinstance(v, Path) else v for v in value]
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> 'ProjectConfig':
        """Create from dictionary, converting strings to Paths."""
        # Remove None values to allow defaults
        d = {k: v for k, v in d.items() if v is not None}
        return cls(**d)


class Project:
    """Represents an IOBRpy analysis project."""

    def __init__(self, name: str, root_dir: Path, config: Optional[ProjectConfig] = None):
        """
        Initialize a project.

        Args:
            name: Project name
            root_dir: Root directory for projects
            config: Optional pre-existing configuration
        """
        self.name = name
        self.root_dir = Path(root_dir)
        self.project_dir = self.root_dir / name
        self.config = config or ProjectConfig(name=name, workspace_dir=self.project_dir)

    @classmethod
    def create(cls, name: str, root_dir: Path, **kwargs) -> 'Project':
        """
        Create a new project.

        Args:
            name: Project name
            root_dir: Root directory for projects
            **kwargs: Additional configuration options

        Returns:
            New Project instance
        """
        project = cls(name, root_dir)
        project.config.workspace_dir = project.project_dir

        # Apply any additional settings
        for key, value in kwargs.items():
            if hasattr(project.config, key):
                setattr(project.config, key, value)
            else:
                project.config.custom_settings[key] = value

        # Create project directory
        project.project_dir.mkdir(parents=True, exist_ok=True)

        # Create standard subdirectories
        (project.project_dir / 'input').mkdir(exist_ok=True)
        (project.project_dir / 'output').mkdir(exist_ok=True)
        (project.project_dir / 'logs').mkdir(exist_ok=True)

        # Save configuration
        project.save_config()

        return project

    @classmethod
    def load(cls, name: str, root_dir: Path) -> 'Project':
        """
        Load an existing project.

        Args:
            name: Project name
            root_dir: Root directory for projects

        Returns:
            Loaded Project instance

        Raises:
            FileNotFoundError: If project doesn't exist
        """
        project_dir = Path(root_dir) / name
        config_file = project_dir / 'project.json'

        if not config_file.exists():
            raise FileNotFoundError(f"Project config not found: {config_file}")

        with open(config_file, 'r') as f:
            config_dict = json.load(f)

        config = ProjectConfig.from_dict(config_dict)
        project = cls(name, root_dir, config=config)
        project.project_dir = project_dir

        return project

    def save_config(self) -> None:
        """Save project configuration to file."""
        config_file = self.project_dir / 'project.json'
        with open(config_file, 'w') as f:
            json.dump(self.config.to_dict(), f, indent=2)

    def update_config(self, **kwargs) -> None:
        """
        Update project configuration.

        Args:
            **kwargs: Configuration options to update
        """
        for key, value in kwargs.items():
            if hasattr(self.config, key):
                setattr(self.config, key, value)
            else:
                self.config.custom_settings[key] = value
        self.save_config()

    def get_input_dir(self) -> Path:
        """Get input directory path."""
        return self.project_dir / 'input'

    def get_output_dir(self) -> Path:
        """Get output directory path."""
        return self.project_dir / 'output'

    def get_logs_dir(self) -> Path:
        """Get logs directory path."""
        return self.project_dir / 'logs'

    def list_inputs(self) -> List[Path]:
        """List all files in input directory."""
        input_dir = self.get_input_dir()
        if not input_dir.exists():
            return []
        return list(input_dir.iterdir())

    def list_outputs(self, pattern: str = '*') -> List[Path]:
        """
        List output files matching pattern.

        Args:
            pattern: Glob pattern for matching files

        Returns:
            List of matching Path objects
        """
        output_dir = self.get_output_dir()
        if not output_dir.exists():
            return []
        return list(output_dir.glob(pattern))

    def get_status(self) -> Dict[str, Any]:
        """
        Get project status information.

        Returns:
            Dictionary with status details
        """
        return {
            'name': self.name,
            'path': str(self.project_dir),
            'created_at': self.config.created_at,
            'input_type': self.config.input_type,
            'organism': self.config.organism,
            'mode': self.config.mode,
            'input_files': len(self.list_inputs()),
            'output_files': len(self.list_outputs()),
        }


class ProjectManager:
    """Manages multiple IOBRpy projects."""

    def __init__(self, root_dir: Path):
        """
        Initialize project manager.

        Args:
            root_dir: Root directory for projects
        """
        self.root_dir = Path(root_dir)
        self.root_dir.mkdir(parents=True, exist_ok=True)

    def create_project(self, name: str, **kwargs) -> Project:
        """
        Create a new project.

        Args:
            name: Project name
            **kwargs: Project configuration options

        Returns:
            New Project instance
        """
        if (self.root_dir / name).exists():
            raise FileExistsError(f"Project '{name}' already exists")

        return Project.create(name, self.root_dir, **kwargs)

    def get_project(self, name: str) -> Project:
        """
        Get an existing project.

        Args:
            name: Project name

        Returns:
            Project instance

        Raises:
            FileNotFoundError: If project doesn't exist
        """
        if not (self.root_dir / name).exists():
            raise FileNotFoundError(f"Project '{name}' not found")
        return Project.load(name, self.root_dir)

    def list_projects(self) -> List[Dict[str, Any]]:
        """
        List all projects.

        Returns:
            List of project info dictionaries
        """
        projects = []
        for path in self.root_dir.iterdir():
            if path.is_dir() and (path / 'project.json').exists():
                try:
                    project = Project.load(path.name, self.root_dir)
                    projects.append(project.get_status())
                except Exception:
                    # Skip invalid projects
                    pass
        return projects

    def delete_project(self, name: str) -> bool:
        """
        Delete a project.

        Args:
            name: Project name

        Returns:
            True if deleted, False if not found

        Raises:
            OSError: If deletion fails for other reasons
        """
        project_dir = self.root_dir / name
        if not project_dir.exists():
            return False

        # Delete project directory
        import shutil
        try:
            shutil.rmtree(project_dir)
            return True
        except OSError:
            return False
