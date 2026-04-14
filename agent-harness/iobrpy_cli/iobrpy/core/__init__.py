"""
core - Core modules for IOBRpy CLI harness.

This package provides the core functionality for managing IOBRpy projects,
sessions, TME analysis, and quantification workflows.
"""

from .project import Project, ProjectConfig, ProjectManager
from .session import Session, SessionContext, SessionManager, format_command_history
from .export import Exporter, ExportResult
from .tme import (
    TMEAnalyzer,
    DeconvolutionMethod,
    ScoringMethod,
    DeconvolutionResult,
    SignatureScoreResult,
    LRResult,
)
from .hla_tcr import (
    HLATCRAnalyzer,
    Trust4Mode,
    Trust4Result,
    HLATypingResult,
    SpecHLAResult,
    ExtractHLAResult,
)
from .quantification import (
    QuantificationWorkflow,
    QuantificationMode,
    QCResult,
    QuantificationResult,
    MergeResult,
    TPMResult,
    CleanResult,
    ValidateResult,
    ScanResult,
    VersionInfo,
)
from .workflow import WorkflowExecutor, WorkflowResult

__all__ = [
    # Project
    'Project',
    'ProjectConfig',
    'ProjectManager',
    # Session
    'Session',
    'SessionContext',
    'SessionManager',
    'format_command_history',
    # Export
    'Exporter',
    'ExportResult',
    # TME
    'TMEAnalyzer',
    'DeconvolutionMethod',
    'ScoringMethod',
    'DeconvolutionResult',
    'SignatureScoreResult',
    'LRResult',
    # HLA/TCR/BCR
    'HLATCRAnalyzer',
    'Trust4Mode',
    'Trust4Result',
    'HLATypingResult',
    'SpecHLAResult',
    'ExtractHLAResult',
    # Quantification
    'QuantificationWorkflow',
    'QuantificationMode',
    'QCResult',
    'QuantificationResult',
    'MergeResult',
    'TPMResult',
    'CleanResult',
    'ValidateResult',
    'ScanResult',
    'VersionInfo',
    # Workflow
    'WorkflowExecutor',
    'WorkflowResult',
]
