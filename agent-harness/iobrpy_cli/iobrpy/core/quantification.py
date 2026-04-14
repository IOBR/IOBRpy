"""
quantification.py - Quantification workflows for IOBRpy CLI harness.

This module provides high-level interfaces for FASTQ-to-TPM quantification
workflows including Salmon and STAR alignment.
"""

import subprocess
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional, List, Dict, Any


class QuantificationMode(Enum):
    """Quantification mode."""

    SALMON = "salmon"
    STAR = "star"


@dataclass
class QCResult:
    """Result of FASTQ QC operation."""

    success: bool
    output_dir: Optional[Path] = None
    sample_count: int = 0
    message: str = ""
    duration_ms: Optional[int] = None


@dataclass
class QuantificationResult:
    """Result of quantification operation."""

    mode: QuantificationMode
    success: bool
    output_dir: Optional[Path] = None
    sample_count: int = 0
    message: str = ""
    duration_ms: Optional[int] = None


@dataclass
class MergeResult:
    """Result of merge operation."""

    success: bool
    output_file: Optional[Path] = None
    message: str = ""


@dataclass
class TPMResult:
    """Result of TPM matrix generation."""

    success: bool
    output_file: Optional[Path] = None
    gene_count: int = 0
    sample_count: int = 0
    message: str = ""


@dataclass
class CleanResult:
    """Result of pipeline clean operation."""

    success: bool
    output_dir: Optional[Path] = None
    files_removed: int = 0
    dirs_removed: int = 0
    space_freed_mb: float = 0.0
    deleted_count: int = 0
    message: str = ""
    duration_ms: Optional[int] = None


@dataclass
class ValidateResult:
    """Result of pipeline validation operation."""

    success: bool
    is_valid: bool = False
    output_dir: Optional[Path] = None
    completed_steps: List[str] = field(default_factory=list)
    missing_steps: List[str] = field(default_factory=list)
    issues: List[str] = field(default_factory=list)
    message: str = ""
    duration_ms: Optional[int] = None


@dataclass
class ScanResult:
    """Result of FASTQ scan operation."""

    success: bool
    sample_count: int = 0
    samples: List[str] = field(default_factory=list)
    issues: List[str] = field(default_factory=list)
    message: str = ""
    duration_ms: Optional[int] = None


@dataclass
class VersionInfo:
    """Version information result."""

    iobrpy_version: str = ""
    iobrpy_installed: bool = False
    external_tools: Optional[Dict[str, str]] = None
    python_version: str = ""
    message: str = ""
    duration_ms: Optional[int] = None


class QuantificationWorkflow:
    """High-level interface for quantification workflows."""

    def __init__(self, iobrpy_cli: str = "iobrpy"):
        """
        Initialize quantification workflow.

        Args:
            iobrpy_cli: Path to iobrpy CLI command
        """
        self.iobrpy_cli = iobrpy_cli

    def run_fastq_qc(
        self,
        fastq_dir: Path,
        output_dir: Path,
        threads: int = 8,
        batch_size: int = 1,
        suffix1: str = "_R1",
        suffix2: Optional[str] = None,
        single_end: bool = False,
        length_required: int = 50,
    ) -> QCResult:
        """
        Run FASTQ quality control with fastp.

        Args:
            fastq_dir: Directory containing raw FASTQ files
            output_dir: Output directory for trimmed FASTQ files
            threads: Number of threads
            batch_size: Batch size for parallel processing
            suffix1: Suffix for R1/fastq files
            suffix2: Suffix for R2 files (paired-end)
            single_end: Single-end reads
            length_required: Minimum read length

        Returns:
            QCResult
        """
        import time

        start_time = time.time()
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            self.iobrpy_cli,
            "fastq_qc",
            "--path1_fastq", str(fastq_dir),
            "--path2_fastp", str(output_dir),
            "--num_threads", str(threads),
            "--batch_size", str(batch_size),
            "--suffix1", suffix1,
            "--length_required", str(length_required),
        ]

        if single_end:
            cmd.append("--se")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )

            duration = int((time.time() - start_time) * 1000)

            if result.returncode == 0:
                task_markers = list(output_dir.glob("*.task.complete"))
                if task_markers:
                    sample_count = len(task_markers)
                elif single_end:
                    sample_count = len([
                        f for f in output_dir.iterdir()
                        if f.is_file() and f.name.endswith(suffix1)
                    ])
                else:
                    inferred_suffix2 = suffix2
                    if inferred_suffix2 is None:
                        if "_R1" in suffix1:
                            inferred_suffix2 = suffix1.replace("_R1", "_R2")
                        elif "_1" in suffix1:
                            inferred_suffix2 = suffix1.replace("_1", "_2", 1)
                        else:
                            inferred_suffix2 = "_R2"
                    sample_count = self.scan_fastq(
                        output_dir,
                        suffix1=suffix1,
                        suffix2=inferred_suffix2,
                        single_end=single_end,
                    ).sample_count

                return QCResult(
                    success=True,
                    output_dir=output_dir,
                    sample_count=sample_count,
                    duration_ms=duration,
                    message=f"QC completed for {sample_count} samples",
                )
            else:
                return QCResult(
                    success=False,
                    message=f"QC failed: {result.stderr}",
                    duration_ms=duration,
                )

        except FileNotFoundError:
            return QCResult(
                success=False,
                message=f"iobrpy CLI not found: {self.iobrpy_cli}",
            )
        except Exception as e:
            return QCResult(
                success=False,
                message=f"QC error: {str(e)}",
            )

    def run_salmon_quantification(
        self,
        fastq_dir: Path,
        output_dir: Path,
        index: Path,
        threads: int = 8,
        batch_size: int = 1,
        suffix1: str = "_R1",
        suffix2: Optional[str] = None,
        gtf: Optional[Path] = None,
    ) -> QuantificationResult:
        """
        Run Salmon quantification.

        Args:
            fastq_dir: Directory containing FASTQ files
            output_dir: Output directory
            index: Salmon index directory
            threads: Number of threads
            batch_size: Batch size for parallel processing
            suffix1: Suffix for R1 files
            suffix2: Suffix for R2 files
            gtf: Optional GTF file for quantification

        Returns:
            QuantificationResult
        """
        import time

        start_time = time.time()
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            self.iobrpy_cli,
            "batch_salmon",
            "--index", str(index),
            "--path_fq", str(fastq_dir),
            "--path_out", str(output_dir),
            "--num_threads", str(threads),
            "--batch_size", str(batch_size),
            "--suffix1", suffix1,
        ]

        if gtf:
            cmd.extend(["--gtf", str(gtf)])

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )

            duration = int((time.time() - start_time) * 1000)

            if result.returncode == 0:
                # Count output directories (one per sample)
                sample_count = len([d for d in output_dir.iterdir() if d.is_dir()])

                return QuantificationResult(
                    mode=QuantificationMode.SALMON,
                    success=True,
                    output_dir=output_dir,
                    sample_count=sample_count,
                    duration_ms=duration,
                    message=f"Salmon quantification completed for {sample_count} samples",
                )
            else:
                return QuantificationResult(
                    mode=QuantificationMode.SALMON,
                    success=False,
                    message=f"Salmon quantification failed: {result.stderr}",
                    duration_ms=duration,
                )

        except FileNotFoundError:
            return QuantificationResult(
                mode=QuantificationMode.SALMON,
                success=False,
                message=f"iobrpy CLI not found: {self.iobrpy_cli}",
            )
        except Exception as e:
            return QuantificationResult(
                mode=QuantificationMode.SALMON,
                success=False,
                message=f"Salmon quantification error: {str(e)}",
            )

    def run_star_alignment(
        self,
        fastq_dir: Path,
        output_dir: Path,
        index: Path,
        threads: int = 8,
        batch_size: int = 1,
        suffix1: str = "_R1",
        suffix2: Optional[str] = None,
    ) -> QuantificationResult:
        """
        Run STAR alignment and counting.

        Args:
            fastq_dir: Directory containing FASTQ files
            output_dir: Output directory
            index: STAR index directory
            threads: Number of threads
            batch_size: Batch size for parallel processing
            suffix1: Suffix for R1 files
            suffix2: Suffix for R2 files

        Returns:
            QuantificationResult
        """
        import time

        start_time = time.time()
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            self.iobrpy_cli,
            "batch_star_count",
            "--index", str(index),
            "--path_fq", str(fastq_dir),
            "--path_out", str(output_dir),
            "--num_threads", str(threads),
            "--batch_size", str(batch_size),
            "--suffix1", suffix1,
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )

            duration = int((time.time() - start_time) * 1000)

            if result.returncode == 0:
                # Count output files
                sample_count = len(list(output_dir.glob("*ReadsPerGene*")))

                return QuantificationResult(
                    mode=QuantificationMode.STAR,
                    success=True,
                    output_dir=output_dir,
                    sample_count=sample_count,
                    duration_ms=duration,
                    message=f"STAR alignment completed for {sample_count} samples",
                )
            else:
                return QuantificationResult(
                    mode=QuantificationMode.STAR,
                    success=False,
                    message=f"STAR alignment failed: {result.stderr}",
                    duration_ms=duration,
                )

        except FileNotFoundError:
            return QuantificationResult(
                mode=QuantificationMode.STAR,
                success=False,
                message=f"iobrpy CLI not found: {self.iobrpy_cli}",
            )
        except Exception as e:
            return QuantificationResult(
                mode=QuantificationMode.STAR,
                success=False,
                message=f"STAR alignment error: {str(e)}",
            )

    def merge_salmon_results(
        self,
        salmon_dir: Path,
        project: str = "merged",
        threads: int = 1,
    ) -> MergeResult:
        """
        Merge Salmon quantification results.

        Args:
            salmon_dir: Directory containing Salmon quantification outputs
            project: Project name for merged output
            threads: Number of threads

        Returns:
            MergeResult
        """
        salmon_dir = Path(salmon_dir)

        cmd = [
            self.iobrpy_cli,
            "merge_salmon",
            "--path_salmon", str(salmon_dir),
            "--project", project,
            "--num_processes", str(threads),
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=str(salmon_dir),
            )

            if result.returncode == 0:
                # Find the merged output file
                for pattern in ["*_salmon_tpm.tsv", "*_salmon_tpm.tsv.gz"]:
                    for f in salmon_dir.glob(pattern):
                        return MergeResult(
                            success=True,
                            output_file=f,
                            message=f"Successfully merged Salmon results to {f.name}",
                        )

                return MergeResult(
                    success=True,
                    message="Salmon merge completed (output file not found in expected location)",
                )
            else:
                return MergeResult(
                    success=False,
                    message=f"Salmon merge failed: {result.stderr}",
                )

        except FileNotFoundError:
            return MergeResult(
                success=False,
                message=f"iobrpy CLI not found: {self.iobrpy_cli}",
            )
        except Exception as e:
            return MergeResult(
                success=False,
                message=f"Salmon merge error: {str(e)}",
            )

    def merge_star_results(
        self,
        star_dir: Path,
        project: str = "merged",
    ) -> MergeResult:
        """
        Merge STAR count results.

        Args:
            star_dir: Directory containing STAR count outputs
            project: Project name for merged output

        Returns:
            MergeResult
        """
        star_dir = Path(star_dir)

        cmd = [
            self.iobrpy_cli,
            "merge_star_count",
            "--path", str(star_dir),
            "--project", project,
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=str(star_dir),
            )

            if result.returncode == 0:
                # Find the merged output file
                for pattern in ["*_star_ReadsPerGene*", "*_ReadsPerGene*"]:
                    for f in star_dir.glob(pattern):
                        return MergeResult(
                            success=True,
                            output_file=f,
                            message=f"Successfully merged STAR results to {f.name}",
                        )

                return MergeResult(
                    success=True,
                    message="STAR merge completed (output file not found in expected location)",
                )
            else:
                return MergeResult(
                    success=False,
                    message=f"STAR merge failed: {result.stderr}",
                )

        except FileNotFoundError:
            return MergeResult(
                success=False,
                message=f"iobrpy CLI not found: {self.iobrpy_cli}",
            )
        except Exception as e:
            return MergeResult(
                success=False,
                message=f"STAR merge error: {str(e)}",
            )

    def prepare_salmon_tpm(
        self,
        input_file: Path,
        output_file: Path,
        return_feature: str = "symbol",
        remove_version: bool = True,
    ) -> TPMResult:
        """
        Prepare Salmon TPM matrix.

        Args:
            input_file: Merged Salmon TPM file
            output_file: Output TPM matrix file
            return_feature: Gene feature type (ENST, ENSG, symbol)
            remove_version: Remove version suffix from gene IDs

        Returns:
            TPMResult
        """
        cmd = [
            self.iobrpy_cli,
            "prepare_salmon",
            "--input", str(input_file),
            "--output", str(output_file),
            "--return_feature", return_feature,
        ]

        if remove_version:
            cmd.append("--remove_version")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )

            if result.returncode == 0 and Path(output_file).exists():
                # Try to read dimensions
                gene_count = 0
                sample_count = 0
                try:
                    import pandas as pd
                    df = pd.read_csv(output_file)
                    gene_count = len(df)
                    sample_count = len(df.columns) - 1 if 'ID' in df.columns else len(df.columns)
                except Exception:
                    pass

                return TPMResult(
                    success=True,
                    output_file=Path(output_file),
                    gene_count=gene_count,
                    sample_count=sample_count,
                    message=f"Prepared TPM matrix: {gene_count} genes, {sample_count} samples",
                )
            else:
                return TPMResult(
                    success=False,
                    message=f"Salmon TPM preparation failed: {result.stderr}",
                )

        except FileNotFoundError:
            return TPMResult(
                success=False,
                message=f"iobrpy CLI not found: {self.iobrpy_cli}",
            )
        except Exception as e:
            return TPMResult(
                success=False,
                message=f"Salmon TPM preparation error: {str(e)}",
            )

    def count2tpm(
        self,
        input_file: Path,
        output_file: Path,
        idtype: str = "ensembl",
        org: str = "hsa",
        source: str = "local",
        remove_version: bool = True,
        **kwargs,
    ) -> TPMResult:
        """
        Convert count matrix to TPM.

        Args:
            input_file: Count matrix file
            output_file: Output TPM matrix file
            idtype: Gene ID type
            org: Organism (hsa, mmus)
            source: Gene length source
            remove_version: Remove version suffix
            **kwargs: Additional parameters

        Returns:
            TPMResult
        """
        cmd = [
            self.iobrpy_cli,
            "count2tpm",
            "--input", str(input_file),
            "--output", str(output_file),
            "--idtype", idtype,
            "--org", org,
            "--source", source,
        ]

        if remove_version:
            cmd.append("--remove_version")

        if kwargs.get("check_data"):
            cmd.append("--check_data")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )

            if result.returncode == 0 and Path(output_file).exists():
                # Try to read dimensions
                gene_count = 0
                sample_count = 0
                try:
                    import pandas as pd
                    df = pd.read_csv(output_file)
                    gene_count = len(df)
                    sample_count = len(df.columns) - 1 if 'ID' in df.columns else len(df.columns)
                except Exception:
                    pass

                return TPMResult(
                    success=True,
                    output_file=Path(output_file),
                    gene_count=gene_count,
                    sample_count=sample_count,
                    message=f"Converted counts to TPM: {gene_count} genes, {sample_count} samples",
                )
            else:
                return TPMResult(
                    success=False,
                    message=f"Count2TPM conversion failed: {result.stderr}",
                )

        except FileNotFoundError:
            return TPMResult(
                success=False,
                message=f"iobrpy CLI not found: {self.iobrpy_cli}",
            )
        except Exception as e:
            return TPMResult(
                success=False,
                message=f"Count2TPM error: {str(e)}",
            )

    def log2_transform(
        self,
        input_file: Path,
        output_file: Path,
    ) -> TPMResult:
        """
        Apply log2(x+1) transformation.

        Args:
            input_file: Input TPM matrix
            output_file: Output log2-transformed matrix

        Returns:
            TPMResult
        """
        cmd = [
            self.iobrpy_cli,
            "log2_eset",
            "-i", str(input_file),
            "-o", str(output_file),
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )

            if result.returncode == 0 and Path(output_file).exists():
                return TPMResult(
                    success=True,
                    output_file=Path(output_file),
                    message="Log2 transformation completed",
                )
            else:
                return TPMResult(
                    success=False,
                    message=f"Log2 transformation failed: {result.stderr}",
                )

        except FileNotFoundError:
            return TPMResult(
                success=False,
                message=f"iobrpy CLI not found: {self.iobrpy_cli}",
            )
        except Exception as e:
            return TPMResult(
                success=False,
                message=f"Log2 transformation error: {str(e)}",
            )

    def run_salmon_pipeline(
        self,
        fastq_dir: Path,
        output_dir: Path,
        salmon_index: Path,
        threads: int = 8,
        batch_size: int = 1,
        suffix1: str = "_R1",
        length_required: int = 50,
        return_feature: str = "symbol",
    ) -> Dict[str, Any]:
        """
        Run complete Salmon quantification pipeline.

        Args:
            fastq_dir: Directory containing raw FASTQ files
            output_dir: Main output directory
            salmon_index: Salmon index path
            threads: Number of threads
            batch_size: Batch size
            suffix1: FASTQ suffix pattern
            length_required: Minimum read length
            return_feature: Gene feature type

        Returns:
            Dictionary with all pipeline results
        """
        output_dir = Path(output_dir)

        results = {
            'qc': None,
            'quantification': None,
            'merge': None,
            'tpm_matrix': None,
            'log2_matrix': None,
        }

        # 1. FASTQ QC
        qc_dir = output_dir / "01-qc"
        results['qc'] = self.run_fastq_qc(
            fastq_dir=fastq_dir,
            output_dir=qc_dir,
            threads=threads,
            batch_size=batch_size,
            suffix1=suffix1,
            length_required=length_required,
        )

        if not results['qc'].success:
            return results

        # 2. Salmon quantification
        salmon_dir = output_dir / "02-salmon"
        results['quantification'] = self.run_salmon_quantification(
            fastq_dir=qc_dir,
            output_dir=salmon_dir,
            index=salmon_index,
            threads=threads,
            batch_size=batch_size,
            suffix1=suffix1,
        )

        if not results['quantification'].success:
            return results

        # 3. Merge results
        results['merge'] = self.merge_salmon_results(
            salmon_dir=salmon_dir,
            project="runall",
            threads=threads,
        )

        if not results['merge'].success:
            return results

        # 4. Prepare TPM matrix
        tpm_dir = output_dir / "03-tpm"
        tpm_dir.mkdir(parents=True, exist_ok=True)

        if results['merge'].output_file:
            results['tpm_matrix'] = self.prepare_salmon_tpm(
                input_file=results['merge'].output_file,
                output_file=tpm_dir / "prepare_salmon.csv",
                return_feature=return_feature,
            )
        else:
            results['tpm_matrix'] = TPMResult(success=False, message="No merged file found")

        if not results['tpm_matrix'].success:
            return results

        # 5. Log2 transform
        results['log2_matrix'] = self.log2_transform(
            input_file=results['tpm_matrix'].output_file,
            output_file=tpm_dir / "tpm_matrix.csv",
        )

        return results

    def run_star_pipeline(
        self,
        fastq_dir: Path,
        output_dir: Path,
        star_index: Path,
        threads: int = 8,
        batch_size: int = 1,
        suffix1: str = "_R1",
        length_required: int = 50,
        idtype: str = "ensembl",
        org: str = "hsa",
    ) -> Dict[str, Any]:
        """
        Run complete STAR quantification pipeline.

        Args:
            fastq_dir: Directory containing raw FASTQ files
            output_dir: Main output directory
            star_index: STAR index path
            threads: Number of threads
            batch_size: Batch size
            suffix1: FASTQ suffix pattern
            length_required: Minimum read length
            idtype: Gene ID type
            org: Organism

        Returns:
            Dictionary with all pipeline results
        """
        output_dir = Path(output_dir)

        results = {
            'qc': None,
            'alignment': None,
            'merge': None,
            'tpm_matrix': None,
            'log2_matrix': None,
        }

        # 1. FASTQ QC
        qc_dir = output_dir / "01-qc"
        results['qc'] = self.run_fastq_qc(
            fastq_dir=fastq_dir,
            output_dir=qc_dir,
            threads=threads,
            batch_size=batch_size,
            suffix1=suffix1,
            length_required=length_required,
        )

        if not results['qc'].success:
            return results

        # 2. STAR alignment
        star_dir = output_dir / "02-star"
        results['alignment'] = self.run_star_alignment(
            fastq_dir=qc_dir,
            output_dir=star_dir,
            index=star_index,
            threads=threads,
            batch_size=batch_size,
            suffix1=suffix1,
        )

        if not results['alignment'].success:
            return results

        # 3. Merge results
        results['merge'] = self.merge_star_results(
            star_dir=star_dir,
            project="runall",
        )

        if not results['merge'].success:
            return results

        # 4. Convert counts to TPM
        tpm_dir = output_dir / "03-tpm"
        tpm_dir.mkdir(parents=True, exist_ok=True)

        if results['merge'].output_file:
            results['tpm_matrix'] = self.count2tpm(
                input_file=results['merge'].output_file,
                output_file=tpm_dir / "count2tpm.csv",
                idtype=idtype,
                org=org,
            )
        else:
            results['tpm_matrix'] = TPMResult(success=False, message="No merged file found")

        if not results['tpm_matrix'].success:
            return results

        # 5. Log2 transform
        results['log2_matrix'] = self.log2_transform(
            input_file=results['tpm_matrix'].output_file,
            output_file=tpm_dir / "tpm_matrix.csv",
        )

        return results

    def runall(
        self,
        fastq_dir: Path,
        output_dir: Path,
        mode: QuantificationMode = QuantificationMode.SALMON,
        index: Optional[Path] = None,
        threads: int = 8,
        batch_size: int = 1,
        suffix1: str = "_R1",
        length_required: int = 50,
        return_feature: str = "symbol",
        idtype: str = "ensembl",
        org: str = "hsa",
        deconvolution_methods: Optional[List[str]] = None,
        include_lr: bool = True,
        signature_groups: str = "all",
        run_trust4: bool = False,
        resume: bool = False,
        dry_run: bool = False,
        **kwargs,
    ) -> Dict[str, Any]:
        """
        Run complete end-to-end pipeline from FASTQ to TME results.

        Args:
            fastq_dir: Directory containing raw FASTQ files
            output_dir: Main output directory
            mode: Quantification mode (salmon or star)
            index: Index path (salmon or star index)
            threads: Number of threads
            batch_size: Batch size
            suffix1: FASTQ suffix pattern
            length_required: Minimum read length
            return_feature: Gene feature type for TPM
            idtype: Gene ID type for count2tpm
            org: Organism
            deconvolution_methods: Deconvolution methods to run
            include_lr: Include ligand-receptor analysis
            signature_groups: Signature groups to score
            run_trust4: Run TRUST4 TCR/BCR analysis
            resume: Resume from existing outputs
            dry_run: Print commands without executing
            **kwargs: Additional parameters

        Returns:
            Dictionary with all pipeline results
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        qc_dir = output_dir / "01-qc"

        results = {
            'steps': {},
            'final': {
                'success': False,
                'message': '',
            },
        }

        # Step markers for resume functionality
        quant_subdir = "02-salmon" if mode == QuantificationMode.SALMON else "02-star"
        done_flags = {
            'fastq_qc': output_dir / "01-qc" / ".fastq_qc.done",
            'quantification': output_dir / quant_subdir / ".quantification.done",
            'merge': output_dir / "03-tpm" / ".merge.done",
            'tme': output_dir / ".tme.done",
        }

        # 1. FASTQ QC
        if resume and done_flags['fastq_qc'].exists():
            results['steps']['fastq_qc'] = {
                'success': True,
                'message': 'Skipped (already completed)',
            }
        elif dry_run:
            results['steps']['fastq_qc'] = {
                'success': True,
                'message': '[dry-run] Would run FASTQ QC',
            }
        else:
            qc_result = self.run_fastq_qc(
                fastq_dir=fastq_dir,
                output_dir=qc_dir,
                threads=threads,
                batch_size=batch_size,
                suffix1=suffix1,
                length_required=length_required,
            )
            results['steps']['fastq_qc'] = {
                'success': qc_result.success,
                'message': qc_result.message,
            }
            if qc_result.success:
                done_flags['fastq_qc'].touch()

        if not results['steps']['fastq_qc']['success'] and not resume:
            results['final']['message'] = 'Pipeline stopped at FASTQ QC step'
            return results

        # 2. Quantification (Salmon or STAR)
        quant_dir = output_dir / ("02-salmon" if mode == QuantificationMode.SALMON else "02-star")
        quant_input_dir = qc_dir if qc_dir.exists() else fastq_dir
        merge_result = None

        # Check index requirement first (applies to both dry_run and actual execution)
        if mode == QuantificationMode.SALMON and not index:
            results['final']['message'] = 'Salmon index path required'
            results['final']['output_dir'] = str(output_dir)
            return results
        elif mode == QuantificationMode.STAR and not index:
            results['final']['message'] = 'STAR index path required'
            results['final']['output_dir'] = str(output_dir)
            return results

        if resume and done_flags['quantification'].exists():
            results['steps']['quantification'] = {
                'success': True,
                'message': 'Skipped (already completed)',
            }
        elif dry_run:
            results['steps']['quantification'] = {
                'success': True,
                'message': f'[dry-run] Would run {mode.value} quantification',
            }
        else:
            if mode == QuantificationMode.SALMON:
                quant_result = self.run_salmon_quantification(
                    fastq_dir=quant_input_dir,
                    output_dir=quant_dir,
                    index=Path(index),
                    threads=threads,
                    batch_size=batch_size,
                    suffix1=suffix1,
                )
            else:  # STAR
                quant_result = self.run_star_alignment(
                    fastq_dir=quant_input_dir,
                    output_dir=quant_dir,
                    index=Path(index),
                    threads=threads,
                    batch_size=batch_size,
                    suffix1=suffix1,
                )
            results['steps']['quantification'] = {
                'success': quant_result.success,
                'message': quant_result.message,
            }
            if quant_result.success:
                done_flags['quantification'].touch()

        if not results['steps']['quantification']['success'] and not resume:
            results['final']['message'] = f'Pipeline stopped at {mode.value} quantification step'
            return results

        # 3. Merge results
        tpm_dir = output_dir / "03-tpm"
        tpm_dir.mkdir(parents=True, exist_ok=True)

        if resume and done_flags['merge'].exists():
            results['steps']['merge'] = {
                'success': True,
                'message': 'Skipped (already completed)',
            }
        elif dry_run:
            results['steps']['merge'] = {
                'success': True,
                'message': f'[dry-run] Would merge {mode.value} results',
            }
        else:
            if mode == QuantificationMode.SALMON:
                merge_result = self.merge_salmon_results(
                    salmon_dir=quant_dir,
                    project="runall",
                    threads=threads,
                )
            else:
                merge_result = self.merge_star_results(
                    star_dir=quant_dir,
                    project="runall",
                )
            results['steps']['merge'] = {
                'success': merge_result.success,
                'message': merge_result.message,
            }
            if merge_result.success:
                done_flags['merge'].touch()

        if not results['steps']['merge']['success'] and not resume:
            results['final']['message'] = 'Pipeline stopped at merge step'
            return results

        # 4. Prepare TPM matrix
        if resume and (tpm_dir / "tpm_matrix.csv").exists():
            results['steps']['tpm_matrix'] = {
                'success': True,
                'message': 'Skipped (already completed)',
            }
        elif dry_run:
            results['steps']['tpm_matrix'] = {
                'success': True,
                'message': '[dry-run] Would prepare TPM matrix',
            }
        else:
            # Find merged output file
            merged_file = None
            if merge_result and merge_result.output_file:
                merged_file = merge_result.output_file

            if mode == QuantificationMode.SALMON:
                for f in quant_dir.glob("*salmon_tpm.tsv*"):
                    if f.is_file() and f.stat().st_size > 0:
                        merged_file = f
                        break
            else:
                for f in quant_dir.glob("*.STAR.count.tsv*"):
                    if f.is_file() and f.stat().st_size > 0:
                        merged_file = f
                        break

            if merged_file:
                if mode == QuantificationMode.SALMON:
                    tpm_result = self.prepare_salmon_tpm(
                        input_file=merged_file,
                        output_file=tpm_dir / "prepare_salmon.csv",
                        return_feature=return_feature,
                    )
                else:
                    tpm_result = self.count2tpm(
                        input_file=merged_file,
                        output_file=tpm_dir / "count2tpm.csv",
                        idtype=idtype,
                        org=org,
                    )
            else:
                tpm_result = TPMResult(success=False, message='No merged file found')

            results['steps']['tpm_matrix'] = {
                'success': tpm_result.success,
                'message': tpm_result.message,
            }

        # Get final TPM matrix for downstream analysis
        final_tpm_file = None
        for f in ["tpm_matrix.csv", "prepare_salmon.csv", "count2tpm.csv"]:
            candidate = tpm_dir / f
            if candidate.exists():
                final_tpm_file = candidate
                break

        results['final']['output_dir'] = str(output_dir)

        if not final_tpm_file and not resume and not dry_run:
            if not resume:
                results['final']['message'] = 'Could not find TPM matrix for downstream analysis'
                return results
            # For resume/dry_run, continue even without TPM file
            final_tpm_file = tpm_dir / "tpm_matrix.csv"  # Placeholder

        if not results['steps']['tpm_matrix']['success'] and not resume:
            results['final']['message'] = 'Pipeline stopped at TPM matrix preparation'
            return results

        # 5. TME analysis (signature scoring, deconvolution, LR)
        if resume and done_flags['tme'].exists():
            results['steps']['tme_analysis'] = {
                'success': True,
                'message': 'Skipped (already completed)',
            }
        elif dry_run:
            results['steps']['tme_analysis'] = {
                'success': True,
                'message': '[dry-run] Would run TME analysis',
            }
        else:
            from .tme import TMEAnalyzer, DeconvolutionMethod
            analyzer = TMEAnalyzer(self.iobrpy_cli)
            parsed_methods = None
            if deconvolution_methods:
                parsed_methods = [
                    TMEAnalyzer.parse_deconvolution_method(method_name)
                    for method_name in deconvolution_methods
                ]

            tme_results = analyzer.run_tme_profile(
                input_matrix=final_tpm_file,
                output_dir=output_dir,
                threads=threads,
                signature_groups=signature_groups,
                deconvolution_methods=parsed_methods,
                include_lr=include_lr,
                lr_cancer_type=kwargs.get('lr_cancer_type', 'pancan'),
                signatures_dirname="04-signatures",
                deconvolution_dirname="05-tme",
                lr_dirname="06-LR_cal",
            )
            signature_success = tme_results['signatures'].success
            deconvolution_success = all(
                result.success for result in tme_results['deconvolution'].values()
            )
            lr_success = (not include_lr) or (tme_results['lr'] is not None and tme_results['lr'].success)
            tme_success = signature_success and deconvolution_success and lr_success

            results['steps']['tme_analysis'] = {
                'success': tme_success,
                'signatures': signature_success,
                'deconvolution': {k: v.success for k, v in tme_results['deconvolution'].items()},
                'lr': tme_results['lr'].success if tme_results['lr'] else False,
                'message': 'TME analysis completed' if tme_success else 'TME analysis completed with failures',
            }

            if tme_success:
                done_flags['tme'].touch()

        if not results['steps']['tme_analysis']['success'] and not dry_run:
            results['final']['message'] = 'Pipeline stopped at TME analysis'
            results['final']['output_dir'] = str(output_dir)
            results['final']['tpm_matrix'] = str(final_tpm_file)
            return results

        # 6. TRUST4 TCR/BCR analysis (optional)
        if run_trust4:
            trust4_dir = output_dir / "07-TCRBCR"
            if resume and any(trust4_dir.glob("*")):
                results['steps']['trust4'] = {
                    'success': True,
                    'message': 'Skipped (already completed)',
                }
            elif dry_run:
                results['steps']['trust4'] = {
                    'success': True,
                    'message': '[dry-run] Would run TRUST4',
                }
            else:
                trust4_dir.mkdir(parents=True, exist_ok=True)
                cmd = [
                    self.iobrpy_cli,
                    "trust4",
                    "--fqdir", str(quant_dir if resume else tpm_dir.parent / "01-qc"),
                    "--out", str(trust4_dir),
                ]
                try:
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                    )
                    results['steps']['trust4'] = {
                        'success': result.returncode == 0,
                        'message': 'TRUST4 completed' if result.returncode == 0 else f'TRUST4 failed: {result.stderr}',
                    }
                except Exception as e:
                    results['steps']['trust4'] = {
                        'success': False,
                        'message': f'TRUST4 error: {str(e)}',
                    }

        required_steps = [
            results['steps']['fastq_qc']['success'],
            results['steps']['quantification']['success'],
            results['steps']['merge']['success'],
            results['steps']['tpm_matrix']['success'],
            results['steps']['tme_analysis']['success'],
        ]
        if run_trust4:
            required_steps.append(results['steps'].get('trust4', {}).get('success', False))

        results['final']['success'] = all(required_steps)
        results['final']['message'] = (
            'Pipeline completed successfully'
            if results['final']['success']
            else 'Pipeline completed with failures'
        )
        results['final']['output_dir'] = str(output_dir)
        results['final']['tpm_matrix'] = str(final_tpm_file)

        return results




    def clean_pipeline(
        self,
        output_dir: Path,
        keep_intermediates: bool = False,
        dry_run: bool = False,
    ) -> CleanResult:
        """
        Clean up intermediate files from pipeline run.

        Args:
            output_dir: Pipeline output directory
            keep_intermediates: Preserve stage directories and remove only
                transient artifacts
            dry_run: If True, show what would be deleted without doing it

        Returns:
            CleanResult with cleaning status
        """
        import time
        import shutil

        start_time = time.time()
        output_dir = Path(output_dir)

        if dry_run:
            return CleanResult(
                success=True,
                output_dir=output_dir,
                deleted_count=0,
                message=f"[dry-run] Would clean {output_dir}",
                duration_ms=int((time.time() - start_time) * 1000),
            )

        files_removed = 0
        dirs_removed = 0
        bytes_removed = 0

        def _safe_unlink(path: Path) -> None:
            nonlocal files_removed, bytes_removed
            try:
                if path.exists():
                    bytes_removed += path.stat().st_size
                    path.unlink()
                    files_removed += 1
            except Exception:
                pass

        def _safe_rmtree(path: Path) -> None:
            nonlocal dirs_removed, bytes_removed
            try:
                if path.exists():
                    bytes_removed += sum(
                        child.stat().st_size
                        for child in path.rglob('*')
                        if child.is_file()
                    )
                    shutil.rmtree(path)
                    dirs_removed += 1
            except Exception:
                pass

        transient_patterns = ['*.tmp', '*.log', '*.done']
        for pattern in transient_patterns:
            for file_path in output_dir.rglob(pattern):
                if file_path.is_file():
                    _safe_unlink(file_path)

        if not keep_intermediates:
            for pattern in ['01-qc', '02-*']:
                for dir_path in output_dir.glob(pattern):
                    if dir_path.is_dir():
                        _safe_rmtree(dir_path)

        return CleanResult(
            success=True,
            output_dir=output_dir,
            files_removed=files_removed,
            dirs_removed=dirs_removed,
            space_freed_mb=round(bytes_removed / (1024 * 1024), 3),
            deleted_count=files_removed + dirs_removed,
            message=(
                f"Cleaned {files_removed} file(s) and {dirs_removed} directorie(s)"
                if not keep_intermediates
                else f"Removed {files_removed} transient file(s)"
            ),
            duration_ms=int((time.time() - start_time) * 1000),
        )

    def validate_pipeline(
        self,
        output_dir: Path,
        mode: QuantificationMode = QuantificationMode.SALMON,
        check_tpm: bool = True,
    ) -> ValidateResult:
        """
        Validate a pipeline output directory for completeness.
        """
        import time

        start_time = time.time()
        output_dir = Path(output_dir)

        completed_steps = []
        missing_steps = []
        issues = []

        expected_dirs = {
            '01-qc': 'FASTQ QC',
            f'02-{mode.value}': f'{mode.value.capitalize()} quantification',
            '03-tpm': 'TPM matrix',
            '04-signatures': 'Signature scores',
            '05-tme': 'Deconvolution',
            '06-LR_cal': 'Ligand-receptor analysis',
        }

        for dir_name, description in expected_dirs.items():
            dir_path = output_dir / dir_name
            if dir_path.exists():
                if any(dir_path.iterdir()):
                    completed_steps.append(description)
                else:
                    issues.append(f"{description} directory exists but is empty")
            else:
                missing_steps.append(description)

        done_flags = {
            '.fastq_qc.done': output_dir / '01-qc' / '.fastq_qc.done',
            '.quantification.done': output_dir / f'02-{mode.value}' / '.quantification.done',
            '.merge.done': output_dir / '03-tpm' / '.merge.done',
            '.tme.done': output_dir / '.tme.done',
        }

        for flag_name, flag_path in done_flags.items():
            if flag_path.exists():
                completed_steps.append(f"Done flag: {flag_name}")

        if check_tpm:
            tpm_dir = output_dir / '03-tpm'
            if tpm_dir.exists():
                tpm_files = []
                for pattern in ('tpm_matrix*.csv', 'prepare_salmon*.csv', 'count2tpm*.csv'):
                    tpm_files.extend(tpm_dir.glob(pattern))
                if tpm_files:
                    completed_steps.append('TPM matrix')
                else:
                    issues.append('TPM matrix directory exists but no matrix found')

        is_valid = len(missing_steps) == 0 and len(issues) == 0
        duration = int((time.time() - start_time) * 1000)

        return ValidateResult(
            success=True,
            is_valid=is_valid,
            output_dir=output_dir,
            completed_steps=completed_steps,
            missing_steps=missing_steps,
            issues=issues,
            message=f"Validation: {'PASSED' if is_valid else 'FAILED'}",
            duration_ms=duration,
        )

    def scan_fastq(
        self,
        fastq_dir: Path,
        suffix1: str = "_R1",
        suffix2: str = "_R2",
        single_end: bool = False,
    ) -> ScanResult:
        """
        Scan FASTQ directory and report structure.
        """
        import time

        start_time = time.time()
        fastq_dir = Path(fastq_dir)

        if not fastq_dir.exists():
            return ScanResult(
                success=False,
                message=f"Directory does not exist: {fastq_dir}",
                duration_ms=int((time.time() - start_time) * 1000),
            )

        samples = []
        issues = []

        def _sample_from_read_suffix(filename: str, suffix: str) -> str | None:
            if filename.endswith(suffix):
                return filename[:-len(suffix)]

            # Users commonly pass read markers such as "_R1" while files keep
            # the FASTQ extension, for example sample_R1.fastq.gz.
            lower_name = filename.lower()
            for ext in (
                ".fastq.gz",
                ".fq.gz",
                ".fastq.bz2",
                ".fq.bz2",
                ".fastq.xz",
                ".fq.xz",
                ".fastq",
                ".fq",
            ):
                if lower_name.endswith(ext):
                    stem = filename[:-len(ext)]
                    if stem.endswith(suffix):
                        return stem[:-len(suffix)]
            return None

        if single_end:
            for f in sorted(fastq_dir.iterdir()):
                if not f.is_file():
                    continue
                sample = _sample_from_read_suffix(f.name, suffix1)
                if sample:
                    samples.append(sample)
        else:
            sample_map = {}

            for f in sorted(fastq_dir.iterdir()):
                if not f.is_file():
                    continue
                sample = _sample_from_read_suffix(f.name, suffix1)
                if sample:
                    sample_map.setdefault(sample, {'r1': None, 'r2': None})['r1'] = f
                    continue
                sample = _sample_from_read_suffix(f.name, suffix2)
                if sample:
                    sample_map.setdefault(sample, {'r1': None, 'r2': None})['r2'] = f

            for sample, files in sample_map.items():
                if files.get('r1') and files.get('r2'):
                    samples.append(sample)
                elif files.get('r1'):
                    issues.append(f"{sample}: Missing R2 file")
                elif files.get('r2'):
                    issues.append(f"{sample}: Missing R1 file")

        duration = int((time.time() - start_time) * 1000)

        return ScanResult(
            success=True,
            sample_count=len(samples),
            samples=samples,
            issues=issues,
            message=f"Found {len(samples)} valid sample pairs" if len(issues) == 0
                     else f"Found {len(samples)} samples with {len(issues)} issues",
            duration_ms=duration,
        )

    def get_version_info(
        self,
        check_external: bool = True,
    ) -> VersionInfo:
        """
        Get version information for IOBRpy and external tools.
        """
        import time
        import sys

        start_time = time.time()

        iobrpy_version = "unknown"
        iobrpy_installed = False
        external_tools = {}

        try:
            result = subprocess.run(
                [self.iobrpy_cli, '--version'],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0:
                iobrpy_installed = True
                iobrpy_version = result.stdout.strip()
        except (FileNotFoundError, subprocess.TimeoutExpired):
            iobrpy_installed = False

        if iobrpy_version == "unknown":
            try:
                import importlib

                iobrpy_module = importlib.import_module("iobrpy")
                module_version = getattr(iobrpy_module, "__version__", "").strip()
                if module_version:
                    iobrpy_version = module_version
            except Exception:
                pass

        if check_external:
            tools_to_check = {
                'fastp': ['fastp', '--version'],
                'salmon': ['salmon', '--version'],
                'star': ['STAR', '--version'],
                'trust4': ['run-trust4', '--help'],
                'fastqc': ['fastqc', '--version'],
            }

            for tool, cmd in tools_to_check.items():
                try:
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                        timeout=5,
                    )
                    if result.returncode == 0:
                        version = result.stdout.split('\n')[0].strip()
                        external_tools[tool] = version
                except (FileNotFoundError, subprocess.TimeoutExpired):
                    external_tools[tool] = "not found"
                except Exception:
                    external_tools[tool] = "error"

        python_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"

        duration = int((time.time() - start_time) * 1000)

        message_parts = []
        message_parts.append(f"IOBRpy: {iobrpy_version} ({'installed' if iobrpy_installed else 'not found'})")
        message_parts.append(f"Python: {python_version}")

        return VersionInfo(
            iobrpy_version=iobrpy_version,
            iobrpy_installed=iobrpy_installed,
            external_tools=external_tools if check_external else {},
            python_version=python_version,
            message=", ".join(message_parts),
            duration_ms=duration,
        )
