import subprocess
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional


class Trust4Mode(Enum):
    SINGLE = "single"
    BATCH = "batch"


@dataclass
class Trust4Result:
    mode: Trust4Mode
    success: bool
    output_dir: Optional[Path] = None
    sample_count: int = 0
    immune_data_file: Optional[Path] = None
    indices_file: Optional[Path] = None
    message: str = ""
    duration_ms: Optional[int] = None


@dataclass
class HLATypingResult:
    success: bool
    output_dir: Optional[Path] = None
    sample_count: int = 0
    hla_result_file: Optional[Path] = None
    message: str = ""
    duration_ms: Optional[int] = None


@dataclass
class SpecHLAResult:
    success: bool
    output_dir: Optional[Path] = None
    sample_count: int = 0
    message: str = ""
    duration_ms: Optional[int] = None


@dataclass
class ExtractHLAResult:
    success: bool
    output_dir: Optional[Path] = None
    sample_id: str = ""
    message: str = ""
    duration_ms: Optional[int] = None


class HLATCRAnalyzer:
    def __init__(self, iobrpy_cli: str = "iobrpy"):
        self.iobrpy_cli = iobrpy_cli

    def run_trust4(self, bam_file=None, bam_dir=None, fastq_dir=None, 
                   output_dir=None, mode=Trust4Mode.SINGLE, threads=8, **kwargs):
        import time
        start_time = time.time()
        
        if not bam_file and not bam_dir and not fastq_dir:
            return Trust4Result(mode=mode, success=False, 
                           message="Must specify BAM file or BAM/FASTQ directory")
        
        output_root = Path(output_dir) if output_dir else None
        if output_root is not None:
            output_root.mkdir(parents=True, exist_ok=True)

        effective_mode = mode
        if bam_dir or fastq_dir:
            effective_mode = Trust4Mode.BATCH

        cmd = [self.iobrpy_cli, "trust4"]
        if effective_mode == Trust4Mode.SINGLE and bam_file:
            cmd.extend(["-b", str(bam_file)])
        elif bam_dir:
            cmd.extend(["-b", str(bam_dir)])
        elif fastq_dir:
            cmd.extend(["--fqdir", str(fastq_dir)])

        if output_root:
            if effective_mode == Trust4Mode.SINGLE:
                sample_source = Path(bam_file) if bam_file else output_root
                prefix = output_root / f"TRUST_{sample_source.stem}"
                cmd.extend(["-o", str(prefix)])
            else:
                cmd.extend(["-o", str(output_root)])
        cmd.extend(["-t", str(threads)])
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            duration = int((time.time() - start_time) * 1000)
            immune_data_file = None
            indices_file = None
            sample_count = 0
            if output_root:
                immune_candidate = output_root / "trust4_immdata.csv"
                indices_candidate = output_root / "trust4_immune_indices.csv"
                if immune_candidate.exists():
                    immune_data_file = immune_candidate
                if indices_candidate.exists():
                    indices_file = indices_candidate
                sample_count = sum(1 for _ in output_root.rglob("*_report.tsv"))
                if result.returncode == 0 and effective_mode == Trust4Mode.SINGLE and sample_count == 0:
                    sample_count = 1
            return Trust4Result(mode=effective_mode, success=result.returncode == 0,
                           output_dir=output_root,
                           sample_count=sample_count,
                           immune_data_file=immune_data_file,
                           indices_file=indices_file,
                           message="TRUST4 completed" if result.returncode == 0 
                                    else f"TRUST4 failed: {result.stderr}",
                           duration_ms=duration)
        except Exception as e:
            return Trust4Result(mode=mode, success=False,
                           message=f"TRUST4 error: {str(e)}")

    def run_hla_typing(self, bam_dir, output_dir, reference="hg38", 
                     threads=8, use_wgs=False, **kwargs):
        import time
        start_time = time.time()
        
        cmd = [self.iobrpy_cli, "hla_typing", "-b", str(bam_dir),
                "-r", reference, "-o", str(output_dir), "-j", str(threads)]
        cmd.extend(["-u", "0" if use_wgs else "1"])
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            duration = int((time.time() - start_time) * 1000)
            output_dir = Path(output_dir)
            spechla_root = output_dir / "SpecHLA"
            sample_count = sum(1 for sample_dir in spechla_root.iterdir() if sample_dir.is_dir()) if spechla_root.exists() else 0
            if result.returncode == 0 and sample_count == 0:
                sample_count = sum(1 for bam_file in Path(bam_dir).glob("*.bam") if bam_file.is_file())
            hla_result_file = output_dir / "hla_result_merged.txt"
            return HLATypingResult(success=result.returncode == 0, output_dir=output_dir,
                                sample_count=sample_count,
                                hla_result_file=hla_result_file if hla_result_file.exists() else None,
                                message="HLA typing completed" if result.returncode == 0 
                                            else f"HLA typing failed: {result.stderr}",
                                duration_ms=duration)
        except Exception as e:
            return HLATypingResult(success=False, message=f"HLA typing error: {str(e)}")

    def run_spechla(self, sample_name: str = "sample", read1: Path | None = None,
                    read2: Path | None = None, output_dir: Path | None = None,
                    threads: int = 8, use_exon: int = 1, **kwargs):
        import time
        start_time = time.time()
        read1 = Path(read1) if read1 is not None else Path("sample_R1.fastq.gz")
        read2 = Path(read2) if read2 is not None else Path("sample_R2.fastq.gz")
        output_dir = Path(output_dir) if output_dir is not None else Path("spechla_output")

        cmd = [
            self.iobrpy_cli, "spechla",
            "-u", str(use_exon),
            "-n", sample_name,
            "-1", str(read1),
            "-2", str(read2),
            "-o", str(output_dir),
            "-j", str(threads),
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            duration = int((time.time() - start_time) * 1000)
            return SpecHLAResult(success=result.returncode == 0,
                                  output_dir=Path(output_dir) if result.returncode == 0 else None,
                                  sample_count=1 if result.returncode == 0 else 0,
                                  message="SpecHLA completed" if result.returncode == 0
                                           else f"SpecHLA failed: {result.stderr}",
                                  duration_ms=duration)
        except Exception as e:
            return SpecHLAResult(success=False, message=f"SpecHLA error: {str(e)}")

    def run_extract_hla_read(self, sample_id: str, bam_path: Path, ref: str,
                           outdir: Path, no_auto_install: bool = False, **kwargs):
        """Extract HLA-related reads from BAM/CRAM files.

        Args:
            sample_id: Sample name (e.g., NA12878)
            bam_path: Path to sorted and indexed BAM/CRAM file
            ref: Reference genome (hg38 or hg19)
            outdir: Output directory for extracted reads
            no_auto_install: Disable automatic installation of missing tools

        Returns:
            ExtractHLAResult with success status and output information
        """
        import time
        import os
        start_time = time.time()

        # Ensure paths are absolute and expanded
        bam_path = Path(os.path.expanduser(bam_path)).resolve()
        outdir = Path(os.path.expanduser(outdir)).resolve()

        # Build command
        cmd = [
            self.iobrpy_cli,
            "extract_hla_read",
            "-s", sample_id,
            "-b", str(bam_path),
            "-r", ref,
            "-o", str(outdir),
        ]

        if no_auto_install:
            cmd.append("--no-auto-install")

        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            duration = int((time.time() - start_time) * 1000)

            if result.returncode == 0:
                message = f"HLA reads extracted for sample '{sample_id}'"
            else:
                message = f"HLA read extraction failed: {result.stderr}"

            return ExtractHLAResult(
                success=result.returncode == 0,
                output_dir=outdir if result.returncode == 0 else None,
                sample_id=sample_id,
                message=message,
                duration_ms=duration,
            )
        except Exception as e:
            return ExtractHLAResult(
                success=False,
                sample_id=sample_id,
                message=f"HLA read extraction error: {str(e)}"
            )
