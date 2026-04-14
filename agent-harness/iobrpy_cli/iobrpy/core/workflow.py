"""
workflow - Individual workflow steps for IOBRpy CLI harness.

This module provides wrappers for individual IOBRpy workflow steps
that are not part of the high-level orchestration commands.
"""

import sys
import tempfile
import io
from contextlib import redirect_stdout
from pathlib import Path
from typing import Optional, Union, List
from dataclasses import dataclass

import pandas as pd

_NATIVE_WORKFLOW_IMPORT_ERROR: Optional[str] = None

try:
    from iobrpy.workflow.prepare_salmon import prepare_salmon_tpm as prepare_salmon_main
    from iobrpy.workflow.count2tpm import count2tpm as count2tpm_main
    from iobrpy.workflow.anno_eset import main as anno_eset_main
    from iobrpy.workflow.calculate_sig_score import calculate_sig_score as calculate_sig_score_main
    from iobrpy.workflow.batch_salmon import main as batch_salmon_main
    from iobrpy.workflow.merge_salmon import main as merge_salmon_main
    from iobrpy.workflow.batch_star_count import main as batch_star_count_main
    from iobrpy.workflow.merge_star_count import main as merge_star_count_main
    from iobrpy.workflow.fastq_qc import main as fastq_qc_main
    from iobrpy.workflow.log2_eset import main as log2_eset_main
    from iobrpy.workflow.tme_cluster import main as tme_cluster_main
    from iobrpy.workflow.nmf import main as nmf_main
    from iobrpy.workflow.mouse2human_eset import main as mouse2human_eset_main
    from iobrpy.workflow.LR_cal import main as LR_cal_main
    from iobrpy.workflow.IPS import main as ips_main
    from iobrpy.workflow.estimate import main as estimate_main
    from iobrpy.workflow.mcpcounter import main as mcpcounter_main
    from iobrpy.utils.print_colorful_message import print_colorful_message
except ImportError as e:
    _NATIVE_WORKFLOW_IMPORT_ERROR = str(e)
    print(f"Warning: Could not import IOBRpy modules: {e}", file=sys.stderr)
    # Create dummy functions for testing
    def prepare_salmon_main(*args, **kwargs): pass
    def count2tpm_main(*args, **kwargs): pass
    def anno_eset_main(*args, **kwargs): pass
    def calculate_sig_score_main(*args, **kwargs): pass
    def batch_salmon_main(*args, **kwargs): pass
    def merge_salmon_main(*args, **kwargs): pass
    def batch_star_count_main(*args, **kwargs): pass
    def merge_star_count_main(*args, **kwargs): pass
    def fastq_qc_main(*args, **kwargs): pass
    def log2_eset_main(*args, **kwargs): pass
    def tme_cluster_main(*args, **kwargs): pass
    def nmf_main(*args, **kwargs): pass
    def mouse2human_eset_main(*args, **kwargs): pass
    def LR_cal_main(*args, **kwargs): pass
    def ips_main(*args, **kwargs): pass
    def estimate_main(*args, **kwargs): pass
    def mcpcounter_main(*args, **kwargs): pass
    def print_colorful_message(*args, **kwargs): pass


@dataclass
class WorkflowResult:
    """Base class for workflow step results."""
    success: bool
    output_path: Optional[Path] = None
    message: str = ""
    duration_ms: Optional[int] = None

    def to_dict(self):
        """Convert to dictionary for JSON output."""
        return {
            'success': self.success,
            'output_path': str(self.output_path) if self.output_path else None,
            'message': self.message,
            'duration_ms': self.duration_ms,
        }


class WorkflowExecutor:
    """Executor for individual workflow steps."""

    def __init__(self, verbose: bool = False):
        self.verbose = verbose

    def _run_with_argv(self, argv: List[str], runner) -> None:
        """Run an argparse-style entrypoint while restoring sys.argv."""
        original_argv = sys.argv[:]
        try:
            sys.argv = argv
            return self._run_callable(runner)
        finally:
            sys.argv = original_argv

    def _run_callable(self, runner):
        """Run a callable while preventing stdout pollution in JSON mode."""
        if self.verbose:
            return runner()

        stdout_buffer = io.StringIO()
        with redirect_stdout(stdout_buffer):
            result = runner()

        buffered_output = stdout_buffer.getvalue()
        if buffered_output:
            print(buffered_output, file=sys.stderr, end="")

        return result

    def _ensure_outputs_exist(self,
                              action: str,
                              paths: List[Path],
                              require_nonempty_dirs: bool = False) -> Optional[WorkflowResult]:
        """Return a failure result when the native workflow produced no tangible output."""
        problems = []

        for path in paths:
            if not path.exists():
                problems.append(f"missing output at {path}")
                continue

            if require_nonempty_dirs and path.is_dir() and not any(path.iterdir()):
                problems.append(f"empty output directory at {path}")

        if not problems:
            return None

        message = f"{action} did not produce expected output: {'; '.join(problems)}"
        if _NATIVE_WORKFLOW_IMPORT_ERROR:
            message += f". Native iobrpy workflow imports failed earlier: {_NATIVE_WORKFLOW_IMPORT_ERROR}"

        return WorkflowResult(success=False, message=message)

    def _materialize_tabular_result(self, result, output_path: Path) -> None:
        """Persist a DataFrame/Series result when a workflow returns data instead of writing files."""
        if isinstance(result, pd.Series):
            result = result.to_frame()

        if not isinstance(result, pd.DataFrame):
            return

        sep = '\t' if output_path.suffix.lower() in ('.txt', '.tsv') else ','
        use_index = not isinstance(result.index, pd.RangeIndex)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(output_path, sep=sep, index=use_index)

    def prepare_salmon(self,
                      input_path: Union[str, Path],
                      output_matrix: Union[str, Path],
                      return_feature: str = 'symbol',
                      remove_version: bool = False) -> WorkflowResult:
        """Prepare Salmon data matrix."""
        input_path = Path(input_path)
        output_matrix = Path(output_matrix)

        if not input_path.exists():
            return WorkflowResult(
                success=False,
                message=f"Input file not found: {input_path}"
            )

        output_matrix.parent.mkdir(parents=True, exist_ok=True)

        try:
            self._run_callable(
                lambda: prepare_salmon_main(
                    eset_path=str(input_path),
                    output_matrix=str(output_matrix),
                    return_feature=return_feature,
                    remove_version=remove_version
                )
            )

            validation_error = self._ensure_outputs_exist(
                "Salmon matrix preparation",
                [output_matrix],
            )
            if validation_error:
                return validation_error

            message = f"Salmon matrix prepared successfully"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")

            return WorkflowResult(
                success=True,
                output_path=output_matrix,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to prepare Salmon matrix: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def count2tpm(self,
                  count_mat: Union[str, Path],
                  output_path: Union[str, Path],
                  effLength_csv: Optional[Union[str, Path]] = None,
                  idtype: str = 'ensembl',
                  org: str = 'hsa',
                  source: str = 'local',
                  id_col: str = 'id',
                  length_col: str = 'eff_length',
                  gene_symbol_col: str = 'symbol',
                  check_data: bool = False,
                  remove_version: bool = False) -> WorkflowResult:
        """Convert count matrix to TPM."""
        count_mat = Path(count_mat)
        output_path = Path(output_path)
        effLength_csv = Path(effLength_csv) if effLength_csv else None

        if not count_mat.exists():
            return WorkflowResult(
                success=False,
                message=f"Count matrix not found: {count_mat}"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Load count matrix
            if count_mat.suffix == '.gz':
                count_df = pd.read_csv(count_mat, sep='\t', index_col=0, compression='gzip')
            else:
                sep = '\t' if count_mat.suffix in ('.tsv', '.tsv.gz') else ','
                count_df = pd.read_csv(count_mat, sep=sep, index_col=0)

            # Load effective length if provided
            eff_df = None
            if effLength_csv and effLength_csv.exists():
                eff_df = pd.read_csv(effLength_csv)

            # Convert to TPM
            tpm_df = self._run_callable(
                lambda: count2tpm_main(
                    count_mat=count_df,
                    anno_grch38=None,
                    anno_gc_vm32=None,
                    idType=idtype,
                    org=org,
                    source=source,
                    remove_version=remove_version,
                    effLength_df=eff_df,
                    id_col=id_col,
                    length_col=length_col,
                    gene_symbol_col=gene_symbol_col,
                    check_data=check_data
                )
            )

            # Save output
            output_path.parent.mkdir(parents=True, exist_ok=True)
            tpm_df.to_csv(output_path)

            message = f"Count matrix converted to TPM at {output_path}"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Shape: {tpm_df.shape}")

            return WorkflowResult(
                success=True,
                output_path=output_path,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to convert count matrix to TPM: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def anno_eset(self,
                  input_path: Union[str, Path],
                  output_path: Union[str, Path],
                  annotation: str,
                  symbol: str = 'symbol',
                  probe: str = 'id',
                  method: str = 'mean',
                  annotation_file: Optional[Union[str, Path]] = None,
                  annotation_key: Optional[str] = None,
                  remove_version: bool = False) -> WorkflowResult:
        """Annotate expression set."""
        input_path = Path(input_path)
        output_path = Path(output_path)
        annotation_file = Path(annotation_file) if annotation_file else None

        if not input_path.exists():
            return WorkflowResult(
                success=False,
                message=f"Input file not found: {input_path}"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Prepare command line arguments
            argv = [
                'anno_eset_main',
                '--input', str(input_path),
                '--output', str(output_path),
                '--annotation', annotation,
                '--symbol', symbol,
                '--probe', probe,
                '--method', method
            ]

            if remove_version:
                argv.append('--remove_version')
            if annotation_file:
                argv.extend(['--annotation-file', str(annotation_file)])
            if annotation_key:
                argv.extend(['--annotation-key', annotation_key])

            self._run_with_argv(argv, anno_eset_main)

            validation_error = self._ensure_outputs_exist(
                "Expression set annotation",
                [output_path],
            )
            if validation_error:
                return validation_error

            message = f"Expression set annotated successfully"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")

            return WorkflowResult(
                success=True,
                output_path=output_path,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to annotate expression set: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def calculate_sig_score(self,
                            input_path: Union[str, Path],
                            output_path: Union[str, Path],
                            sig_group: List[str],
                            score_method: str = 'pca',
                            mini_gene_count: int = 3,
                            adjust_eset: bool = False,
                            parallel_size: int = 1) -> WorkflowResult:
        """Calculate signature scores."""
        input_path = Path(input_path)
        output_path = Path(output_path)

        if not input_path.exists():
            return WorkflowResult(
                success=False,
                message=f"Input file not found: {input_path}"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Load input matrix
            ext = input_path.suffix.lower()
            if ext == '.csv':
                eset_df = pd.read_csv(input_path, sep=',', index_col=0)
            elif ext == '.txt':
                eset_df = pd.read_csv(input_path, sep='\t', index_col=0)
            else:
                eset_df = pd.read_csv(input_path, sep=None, engine='python', index_col=0)

            # Calculate scores
            scores_df = self._run_callable(
                lambda: calculate_sig_score_main(
                    eset_df,
                    sig_group,
                    score_method,
                    mini_gene_count,
                    adjust_eset,
                    parallel_size
                )
            )

            # Save output
            scores_df.to_csv(output_path, index=False)

            message = f"Signature scores calculated and saved to {output_path}"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Shape: {scores_df.shape}")
                print(f"  Signatures: {len(sig_group)} groups")

            return WorkflowResult(
                success=True,
                output_path=output_path,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to calculate signature scores: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def batch_salmon(self,
                     index: Union[str, Path],
                     path_fq: Union[str, Path],
                     path_out: Union[str, Path],
                     suffix1: str = '_1.fastq.gz',
                     batch_size: int = 1,
                     num_threads: int = 8,
                     gtf: Optional[Union[str, Path]] = None) -> WorkflowResult:
        """Batch-run Salmon quantification."""
        index = Path(index)
        path_fq = Path(path_fq)
        path_out = Path(path_out)
        gtf = Path(gtf) if gtf else None

        if not index.exists():
            return WorkflowResult(
                success=False,
                message=f"Salmon index not found: {index}"
            )

        if not path_fq.exists():
            return WorkflowResult(
                success=False,
                message=f"FASTQ directory not found: {path_fq}"
            )

        path_out.mkdir(parents=True, exist_ok=True)

        try:
            # Prepare command line arguments
            argv = [
                'batch_salmon_main',
                '--index', str(index),
                '--path_fq', str(path_fq),
                '--path_out', str(path_out),
                '--suffix1', suffix1,
                '--batch_size', str(batch_size),
                '--num_threads', str(num_threads)
            ]

            if gtf:
                argv.extend(['--gtf', str(gtf)])

            self._run_with_argv(argv, batch_salmon_main)

            validation_error = self._ensure_outputs_exist(
                "Batch Salmon quantification",
                [path_out],
                require_nonempty_dirs=True,
            )
            if validation_error:
                return validation_error

            message = f"Batch Salmon quantification completed"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Output directory: {path_out}")

            return WorkflowResult(
                success=True,
                output_path=path_out,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to run batch Salmon: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def merge_salmon(self,
                    path_salmon: Union[str, Path],
                    project: str,
                    num_processes: Optional[int] = None) -> WorkflowResult:
        """Merge Salmon quant.sf files."""
        path_salmon = Path(path_salmon)

        if not path_salmon.exists():
            return WorkflowResult(
                success=False,
                message=f"Salmon directory not found: {path_salmon}"
            )

        try:
            # Prepare command line arguments
            argv = [
                'merge_salmon_main',
                '--path_salmon', str(path_salmon),
                '--project', project
            ]

            if num_processes:
                argv.extend(['--num_processes', str(num_processes)])

            self._run_with_argv(argv, merge_salmon_main)

            # The native workflow writes gzipped TSV matrices.
            tpm_file = path_salmon / f"{project}_salmon_tpm.tsv.gz"
            reads_file = path_salmon / f"{project}_salmon_count.tsv.gz"

            validation_error = self._ensure_outputs_exist(
                "Salmon result merge",
                [tpm_file, reads_file],
            )
            if validation_error:
                return validation_error

            message = f"Salmon results merged successfully"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  TPM matrix: {tpm_file}")
                print(f"  NumReads matrix: {reads_file}")

            return WorkflowResult(
                success=True,
                output_path=path_salmon,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to merge Salmon results: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def batch_star_count(self,
                        index: Union[str, Path],
                        path_fq: Union[str, Path],
                        path_out: Union[str, Path],
                        suffix1: str = '_1.fastq.gz',
                        batch_size: int = 1,
                        num_threads: int = 8) -> WorkflowResult:
        """Batch-run STAR with GeneCounts."""
        index = Path(index)
        path_fq = Path(path_fq)
        path_out = Path(path_out)

        if not index.exists():
            return WorkflowResult(
                success=False,
                message=f"STAR index not found: {index}"
            )

        if not path_fq.exists():
            return WorkflowResult(
                success=False,
                message=f"FASTQ directory not found: {path_fq}"
            )

        path_out.mkdir(parents=True, exist_ok=True)

        try:
            # Prepare command line arguments
            argv = [
                'batch_star_count_main',
                '--index', str(index),
                '--path_fq', str(path_fq),
                '--path_out', str(path_out),
                '--suffix1', suffix1,
                '--batch_size', str(batch_size),
                '--num_threads', str(num_threads)
            ]

            self._run_with_argv(argv, batch_star_count_main)

            validation_error = self._ensure_outputs_exist(
                "Batch STAR count",
                [path_out],
                require_nonempty_dirs=True,
            )
            if validation_error:
                return validation_error

            message = f"Batch STAR count completed"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Output directory: {path_out}")

            return WorkflowResult(
                success=True,
                output_path=path_out,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to run batch STAR count: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def merge_star_count(self,
                        path: Union[str, Path],
                        project: str) -> WorkflowResult:
        """Merge STAR *_ReadsPerGene.out.tab files."""
        path = Path(path)

        if not path.exists():
            return WorkflowResult(
                success=False,
                message=f"STAR output directory not found: {path}"
            )

        try:
            # Prepare command line arguments
            argv = [
                'merge_star_count_main',
                '--path', str(path),
                '--project', project
            ]

            self._run_with_argv(argv, merge_star_count_main)

            # The native workflow writes a gzipped TSV matrix.
            output_file = path / f"{project}.STAR.count.tsv.gz"

            validation_error = self._ensure_outputs_exist(
                "STAR count merge",
                [output_file],
            )
            if validation_error:
                return validation_error

            message = f"STAR count results merged successfully"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Output file: {output_file}")

            return WorkflowResult(
                success=True,
                output_path=output_file,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to merge STAR count results: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def fastq_qc(self,
                  path1_fastq: Union[str, Path],
                  path2_fastp: Union[str, Path],
                  num_threads: int = 8,
                  suffix1: str = '_1.fastq.gz',
                  batch_size: int = 1,
                  se: bool = False,
                  length_required: int = 50) -> WorkflowResult:
        """Run FASTQ quality control."""
        path1_fastq = Path(path1_fastq)
        path2_fastp = Path(path2_fastp)

        if not path1_fastq.exists():
            return WorkflowResult(
                success=False,
                message=f"FASTQ directory not found: {path1_fastq}"
            )

        path2_fastp.mkdir(parents=True, exist_ok=True)

        try:
            # Prepare command line arguments
            argv = [
                'fastq_qc_main',
                '--path1_fastq', str(path1_fastq),
                '--path2_fastp', str(path2_fastp),
                '--num_threads', str(num_threads),
                '--suffix1', suffix1,
                '--batch_size', str(batch_size),
                '--length_required', str(length_required)
            ]

            if se:
                argv.append('--se')

            self._run_with_argv(argv, fastq_qc_main)

            validation_error = self._ensure_outputs_exist(
                "FASTQ QC",
                [path2_fastp],
                require_nonempty_dirs=True,
            )
            if validation_error:
                return validation_error

            message = f"FASTQ QC completed"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Output directory: {path2_fastp}")

            return WorkflowResult(
                success=True,
                output_path=path2_fastp,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to run FASTQ QC: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def log2_eset(self,
                  input_path: Union[str, Path],
                  output_path: Union[str, Path]) -> WorkflowResult:
        """Apply log2(x+1) to expression matrix."""
        input_path = Path(input_path)
        output_path = Path(output_path)

        if not input_path.exists():
            return WorkflowResult(
                success=False,
                message=f"Input file not found: {input_path}"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Prepare command line arguments
            argv = [
                'log2_eset_main',
                '-i', str(input_path),
                '-o', str(output_path)
            ]

            self._run_with_argv(argv, log2_eset_main)

            validation_error = self._ensure_outputs_exist(
                "Log2 transformation",
                [output_path],
            )
            if validation_error:
                return validation_error

            message = f"Log2 transformation completed"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Output file: {output_path}")

            return WorkflowResult(
                success=True,
                output_path=output_path,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to apply log2 transformation: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def tme_cluster(self,
                    input_path: Union[str, Path],
                    output_path: Union[str, Path],
                    features: Optional[str] = None,
                    pattern: Optional[str] = None,
                    id_col: Optional[str] = None,
                    scale: bool = True,
                    min_nc: int = 2,
                    max_nc: int = 6,
                    max_iter: int = 10,
                    tol: float = 1e-4,
                    print_intermediate: bool = False,
                    input_sep: Optional[str] = None,
                    output_sep: Optional[str] = None) -> WorkflowResult:
        """Run TME clustering."""
        input_path = Path(input_path)
        output_path = Path(output_path)

        if not input_path.exists():
            return WorkflowResult(
                success=False,
                message=f"Input file not found: {input_path}"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Prepare command line arguments
            argv = [
                'tme_cluster_main',
                '--input', str(input_path),
                '--output', str(output_path)
            ]

            if features:
                argv.extend(['--features', features])
            if pattern:
                argv.extend(['--pattern', pattern])
            if id_col:
                argv.extend(['--id', id_col])

            if scale:
                argv.append('--scale')
            else:
                argv.append('--no-scale')

            argv.extend(['--min_nc', str(min_nc)])
            argv.extend(['--max_nc', str(max_nc)])
            argv.extend(['--max_iter', str(max_iter)])
            argv.extend(['--tol', str(tol)])

            if print_intermediate:
                argv.append('--print_result')

            if input_sep:
                argv.extend(['--input_sep', input_sep])
            if output_sep:
                argv.extend(['--output_sep', output_sep])

            self._run_with_argv(argv, tme_cluster_main)

            validation_error = self._ensure_outputs_exist(
                "TME clustering",
                [output_path],
            )
            if validation_error:
                return validation_error

            message = f"TME clustering completed"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Output file: {output_path}")

            return WorkflowResult(
                success=True,
                output_path=output_path,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to run TME clustering: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def nmf_cluster(self,
                    input_path: Union[str, Path],
                    output_path: Union[str, Path],
                    kmin: int = 2,
                    kmax: int = 8,
                    features: Optional[str] = None,
                    log1p: bool = False,
                    normalize: bool = False,
                    shift: Optional[float] = None,
                    random_state: int = 42,
                    max_iter: int = 1000,
                    skip_k_2: bool = False) -> WorkflowResult:
        """Run NMF-based clustering."""
        input_path = Path(input_path)
        output_path = Path(output_path)

        if not input_path.exists():
            return WorkflowResult(
                success=False,
                message=f"Input file not found: {input_path}"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Prepare command line arguments
            argv = [
                'nmf_main',
                '--input', str(input_path),
                '--output', str(output_path)
            ]

            argv.extend(['--kmin', str(kmin)])
            argv.extend(['--kmax', str(kmax)])

            if features:
                argv.extend(['--features', features])

            if log1p:
                argv.append('--log1p')

            if normalize:
                argv.append('--normalize')

            if shift is not None:
                argv.extend(['--shift', str(shift)])

            argv.extend(['--random-state', str(random_state)])
            argv.extend(['--max-iter', str(max_iter)])

            if skip_k_2:
                argv.append('--skip_k_2')

            self._run_with_argv(argv, nmf_main)

            validation_error = self._ensure_outputs_exist(
                "NMF clustering",
                [output_path],
                require_nonempty_dirs=True,
            )
            if validation_error:
                return validation_error

            message = f"NMF clustering completed"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Output directory: {output_path}")

            return WorkflowResult(
                success=True,
                output_path=output_path,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to run NMF clustering: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def mouse2human(self,
                    input_path: Union[str, Path],
                    output_path: Union[str, Path],
                    is_matrix: bool = False,
                    column_of_symbol: Optional[str] = None,
                    verbose: bool = False,
                    sep: Optional[str] = None,
                    out_sep: Optional[str] = None,
                    progress: bool = True) -> WorkflowResult:
        """Convert mouse gene symbols to human symbols."""
        input_path = Path(input_path)
        output_path = Path(output_path)

        if not input_path.exists():
            return WorkflowResult(
                success=False,
                message=f"Input file not found: {input_path}"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Prepare command line arguments
            argv = [
                'mouse2human_eset_main',
                '-i', str(input_path),
                '-o', str(output_path)
            ]

            if is_matrix:
                argv.append('--is_matrix')

            if column_of_symbol:
                argv.extend(['--column_of_symbol', column_of_symbol])

            if verbose:
                argv.append('--verbose')

            if sep:
                argv.extend(['--sep', sep])

            if out_sep:
                argv.extend(['--out_sep', out_sep])

            if progress:
                argv.append('--progress')

            self._run_with_argv(argv, mouse2human_eset_main)

            validation_error = self._ensure_outputs_exist(
                "Mouse-to-human conversion",
                [output_path],
            )
            if validation_error:
                return validation_error

            message = f"Mouse to human conversion completed"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Output file: {output_path}")

            return WorkflowResult(
                success=True,
                output_path=output_path,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to convert mouse to human symbols: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def lr_cal(self,
               input_path: Union[str, Path],
               output_path: Union[str, Path],
               data_type: str = 'tpm',
               id_type: str = 'ensembl',
               cancer_type: str = 'pancan',
               verbose: bool = False) -> WorkflowResult:
        """Calculate ligand-receptor interactions."""
        input_path = Path(input_path)
        output_path = Path(output_path)

        if not input_path.exists():
            return WorkflowResult(
                success=False,
                message=f"Input file not found: {input_path}"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Prepare command line arguments
            argv = [
                'LR_cal_main',
                '--input', str(input_path),
                '--output', str(output_path),
                '--data_type', data_type,
                '--id_type', id_type,
                '--cancer_type', cancer_type
            ]

            if verbose:
                argv.append('--verbose')

            self._run_with_argv(argv, LR_cal_main)

            validation_error = self._ensure_outputs_exist(
                "Ligand-receptor analysis",
                [output_path],
            )
            if validation_error:
                return validation_error

            message = f"Ligand-receptor analysis completed"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Output file: {output_path}")

            return WorkflowResult(
                success=True,
                output_path=output_path,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to calculate ligand-receptor interactions: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def remove_version(self,
                      input_path: Union[str, Path],
                      output_path: Union[str, Path],
                      is_matrix: bool = False,
                      column_of_symbol: Optional[str] = None,
                      verbose: bool = False,
                      sep: Optional[str] = None,
                      out_sep: Optional[str] = None) -> WorkflowResult:
        """Remove version suffixes from gene IDs."""
        from iobrpy.utils.remove_version import strip_versions_in_dataframe

        input_path = Path(input_path)
        output_path = Path(output_path)

        if not input_path.exists():
            return WorkflowResult(
                success=False,
                message=f"Input file not found: {input_path}"
            )

        if not is_matrix and not column_of_symbol:
            return WorkflowResult(
                success=False,
                message="column_of_symbol is required when is_matrix is False"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            read_kwargs = {}
            if sep is not None:
                read_kwargs["sep"] = sep
            else:
                read_kwargs["sep"] = None
                read_kwargs["engine"] = "python"

            if is_matrix:
                read_kwargs["index_col"] = 0

            df = pd.read_csv(input_path, **read_kwargs)

            # Strip versions
            df_stripped, n_stripped = strip_versions_in_dataframe(
                df,
                column=column_of_symbol if not is_matrix else None,
                on_index=is_matrix,
            )

            # Handle duplicates
            if n_stripped > 0:
                from iobrpy.utils.remove_version import deduplicate_after_stripping
                df_stripped = deduplicate_after_stripping(
                    df_stripped,
                    how='mean',
                    on_index=is_matrix,
                    column=column_of_symbol if not is_matrix else None,
                )

            write_sep = out_sep if out_sep is not None else (
                '\t' if output_path.suffix.lower() in ('.tsv', '.txt') else ','
            )

            if is_matrix:
                df_to_write = df_stripped
                write_index = True
            else:
                index_name = df_stripped.index.name or column_of_symbol or 'symbol'
                df_to_write = df_stripped.reset_index().rename(columns={index_name: column_of_symbol})
                write_index = False

            # Write output
            df_to_write.to_csv(output_path, sep=write_sep, index=write_index)

            message = f"Removed version suffixes from {n_stripped} gene IDs"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Output file: {output_path}")

            return WorkflowResult(
                success=True,
                output_path=output_path,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to remove version suffixes: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def mouse2human_eset(self,
                         input_path: Union[str, Path],
                         output_path: Union[str, Path],
                         is_matrix: bool = False,
                         column_of_symbol: Optional[str] = None,
                         verbose: bool = False,
                         sep: Optional[str] = None,
                         out_sep: Optional[str] = None,
                         progress: bool = True) -> WorkflowResult:
        """Convert mouse gene symbols to human symbols for expression set."""
        input_path = Path(input_path)
        output_path = Path(output_path)

        if not input_path.exists():
            return WorkflowResult(
                success=False,
                message=f"Input file not found: {input_path}"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Prepare command line arguments
            argv = [
                'mouse2human_eset_main',
                '--input', str(input_path),
                '--output', str(output_path),
            ]

            if is_matrix:
                argv.append('--is_matrix')
            if column_of_symbol:
                argv.extend(['--column_of_symbol', column_of_symbol])
            if verbose:
                argv.append('--verbose')
            if sep:
                argv.extend(['--sep', sep])
            if out_sep:
                argv.extend(['--out_sep', out_sep])
            if progress:
                argv.append('--progress')

            self._run_with_argv(argv, mouse2human_eset_main)

            validation_error = self._ensure_outputs_exist(
                "Mouse-to-human conversion",
                [output_path],
            )
            if validation_error:
                return validation_error

            message = f"Mouse to human conversion completed"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Output file: {output_path}")

            return WorkflowResult(
                success=True,
                output_path=output_path,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to convert mouse to human symbols: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)
    def run_ips(self,
               input_path: Union[str, Path],
               output_path: Union[str, Path]) -> WorkflowResult:
        """Calculate Immunophenoscore (IPS)."""
        input_path = Path(input_path)
        output_path = Path(output_path)

        if not input_path.exists():
            return WorkflowResult(
                success=False,
                message=f"Input file not found: {input_path}"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Prepare command line arguments
            argv = [
                'ips_main',
                '--input', str(input_path),
                '--output', str(output_path)
            ]

            result_df = self._run_with_argv(argv, ips_main)
            if not output_path.exists():
                self._materialize_tabular_result(result_df, output_path)

            validation_error = self._ensure_outputs_exist(
                "IPS calculation",
                [output_path],
            )
            if validation_error:
                return validation_error

            message = f"IPS calculation completed"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Output file: {output_path}")

            return WorkflowResult(
                success=True,
                output_path=output_path,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to calculate IPS: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def run_estimate(self,
                     input_path: Union[str, Path],
                     output_path: Union[str, Path],
                     platform: str = 'affymetrix') -> WorkflowResult:
        """Calculate ESTIMATE score."""
        input_path = Path(input_path)
        output_path = Path(output_path)

        if not input_path.exists():
            return WorkflowResult(
                success=False,
                message=f"Input file not found: {input_path}"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Prepare command line arguments
            argv = [
                'estimate_main',
                '--input', str(input_path),
                '--output', str(output_path),
                '--platform', platform
            ]

            result_df = self._run_with_argv(argv, estimate_main)
            if not output_path.exists():
                self._materialize_tabular_result(result_df, output_path)

            validation_error = self._ensure_outputs_exist(
                "ESTIMATE calculation",
                [output_path],
            )
            if validation_error:
                return validation_error

            message = f"ESTIMATE calculation completed"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Platform: {platform}")
                print(f"  Output file: {output_path}")

            return WorkflowResult(
                success=True,
                output_path=output_path,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to calculate ESTIMATE: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)

    def run_mcpcounter(self,
                      input_path: Union[str, Path],
                      output_path: Union[str, Path],
                      features: str) -> WorkflowResult:
        """Run MCPcounter estimation."""
        input_path = Path(input_path)
        output_path = Path(output_path)

        if not input_path.exists():
            return WorkflowResult(
                success=False,
                message=f"Input file not found: {input_path}"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Prepare command line arguments
            argv = [
                'mcpcounter_main',
                '--input', str(input_path),
                '--output', str(output_path),
                '--features', features
            ]

            result_df = self._run_with_argv(argv, mcpcounter_main)
            if not output_path.exists():
                self._materialize_tabular_result(result_df, output_path)

            validation_error = self._ensure_outputs_exist(
                "MCPcounter estimation",
                [output_path],
            )
            if validation_error:
                return validation_error

            message = f"MCPcounter estimation completed"
            if self.verbose:
                print_colorful_message(f"✓ {message}", "green")
                print(f"  Features: {features}")
                print(f"  Output file: {output_path}")

            return WorkflowResult(
                success=True,
                output_path=output_path,
                message=message
            )

        except Exception as e:
            error_msg = f"Failed to run MCPcounter: {str(e)}"
            if self.verbose:
                print_colorful_message(f"✗ {error_msg}", "red")
            return WorkflowResult(success=False, message=error_msg)
