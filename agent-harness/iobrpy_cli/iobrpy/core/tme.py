"""
tme.py - Tumor Microenvironment analysis for IOBRpy CLI harness.

This module provides high-level interfaces for TME analysis operations
including deconvolution, signature scoring, and LR analysis.
"""

import subprocess
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional, List, Dict, Any, Union


class DeconvolutionMethod(Enum):
    """Supported deconvolution methods."""

    CIBERSORT = "cibersort"
    IPS = "IPS"
    ESTIMATE = "estimate"
    MCPCOUNTER = "mcpcounter"
    QUANTISEQ = "quantiseq"
    EPIC = "epic"
    BAYESPRISM = "bayesprism"


class ScoringMethod(Enum):
    """Supported signature scoring methods."""

    PCA = "pca"
    ZSCORE = "zscore"
    SSGSEA = "ssgsea"
    INTEGRATION = "integration"


@dataclass
class DeconvolutionResult:
    """Result of a deconvolution operation."""

    method: DeconvolutionMethod
    success: bool
    output_file: Optional[Path] = None
    message: str = ""
    duration_ms: Optional[int] = None


@dataclass
class SignatureScoreResult:
    """Result of signature scoring."""

    method: ScoringMethod
    success: bool
    output_file: Optional[Path] = None
    signature_count: int = 0
    message: str = ""


@dataclass
class LRResult:
    """Result of ligand-receptor analysis."""

    success: bool
    output_file: Optional[Path] = None
    pair_count: int = 0
    message: str = ""


class TMEAnalyzer:
    """High-level interface for TME analysis."""

    def __init__(self, iobrpy_cli: str = "iobrpy"):
        """
        Initialize TME analyzer.

        Args:
            iobrpy_cli: Path to iobrpy CLI command
        """
        self.iobrpy_cli = iobrpy_cli

    @staticmethod
    def parse_deconvolution_method(name: str) -> DeconvolutionMethod:
        """Normalize a user-facing method name to a supported enum value."""
        normalized = name.strip().lower()
        method_map = {
            "cibersort": DeconvolutionMethod.CIBERSORT,
            "ips": DeconvolutionMethod.IPS,
            "estimate": DeconvolutionMethod.ESTIMATE,
            "mcpcounter": DeconvolutionMethod.MCPCOUNTER,
            "quantiseq": DeconvolutionMethod.QUANTISEQ,
            "epic": DeconvolutionMethod.EPIC,
            "bayesprism": DeconvolutionMethod.BAYESPRISM,
        }
        if normalized not in method_map:
            raise ValueError(f"Unsupported deconvolution method: {name}")
        return method_map[normalized]

    @staticmethod
    def parse_signature_groups(signature_groups: Union[str, List[str]]) -> List[str]:
        """Normalize signature input into a flat list of group names."""
        if isinstance(signature_groups, list):
            groups = signature_groups
        else:
            groups = signature_groups.replace(",", " ").split()
        return [group.strip() for group in groups if str(group).strip()]

    def run_deconvolution(
        self,
        input_matrix: Path,
        output_file: Path,
        method: DeconvolutionMethod = DeconvolutionMethod.CIBERSORT,
        threads: int = 1,
        **kwargs,
    ) -> DeconvolutionResult:
        """
        Run a deconvolution method.

        Args:
            input_matrix: Input expression matrix (TPM)
            output_file: Output file path
            method: Deconvolution method
            threads: Number of threads
            **kwargs: Method-specific parameters

        Returns:
            DeconvolutionResult
        """
        import time

        start_time = time.time()
        kwargs = dict(kwargs)

        # Accept both CLI-style and internal key names.
        if "QN" in kwargs and "qn" not in kwargs:
            kwargs["qn"] = kwargs.pop("QN")

        cmd = [self.iobrpy_cli, method.value, "--input", str(input_matrix), "--output", str(output_file)]

        # Add method-specific parameters
        if method == DeconvolutionMethod.CIBERSORT:
            cmd.extend(["--threads", str(threads)])
            if "perm" in kwargs:
                cmd.extend(["--perm", str(kwargs["perm"])])
            if "qn" in kwargs:
                cmd.extend(["--QN", str(kwargs["qn"]).lower()])
            if "absolute" in kwargs:
                cmd.extend(["--absolute", str(kwargs["absolute"]).lower()])
            if "abs_method" in kwargs:
                cmd.extend(["--abs_method", kwargs["abs_method"]])

        elif method == DeconvolutionMethod.IPS:
            pass  # IPS has no additional parameters

        elif method == DeconvolutionMethod.ESTIMATE:
            cmd.extend(["--platform", kwargs.get("platform", "affymetrix")])

        elif method == DeconvolutionMethod.MCPCOUNTER:
            cmd.extend(["--features", kwargs.get("features", "HUGO_symbols")])

        elif method == DeconvolutionMethod.QUANTISEQ:
            if kwargs.get("arrays"):
                cmd.append("--arrays")
            if kwargs.get("tumor"):
                cmd.append("--tumor")
            if kwargs.get("scale_mrna"):
                cmd.append("--scale_mrna")
            if "method" in kwargs:
                cmd.extend(["--method", kwargs["method"]])
            if "signame" in kwargs:
                cmd.extend(["--signame", kwargs["signame"]])
            if "rmgenes" in kwargs:
                cmd.extend(["--rmgenes", kwargs["rmgenes"]])

        elif method == DeconvolutionMethod.EPIC:
            cmd.extend(["--reference", kwargs.get("reference", "TRef")])

        elif method == DeconvolutionMethod.BAYESPRISM:
            if "sc_dat" in kwargs:
                cmd.extend(["--sc_dat", str(kwargs["sc_dat"])])
            if "cell_state_labels" in kwargs:
                cmd.extend(["--cell_state_labels", str(kwargs["cell_state_labels"])])
            if "cell_type_labels" in kwargs:
                cmd.extend(["--cell_type_labels", str(kwargs["cell_type_labels"])])
            if "key" in kwargs:
                cmd.extend(["--key", kwargs["key"]])

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )

            duration = int((time.time() - start_time) * 1000)

            if result.returncode == 0:
                output_file = Path(output_file)
                return DeconvolutionResult(
                    method=method,
                    success=True,
                    output_file=output_file if output_file.exists() else None,
                    duration_ms=duration,
                    message=f"{method.value} completed successfully",
                )
            else:
                return DeconvolutionResult(
                    method=method,
                    success=False,
                    message=f"{method.value} failed: {result.stderr}",
                    duration_ms=duration,
                )

        except FileNotFoundError:
            return DeconvolutionResult(
                method=method,
                success=False,
                message=f"iobrpy CLI not found: {self.iobrpy_cli}",
            )
        except Exception as e:
            return DeconvolutionResult(
                method=method,
                success=False,
                message=f"{method.value} error: {str(e)}",
            )

    def run_all_deconvolution(
        self,
        input_matrix: Path,
        output_dir: Path,
        methods: Optional[List[DeconvolutionMethod]] = None,
        threads: int = 1,
        **kwargs,
    ) -> Dict[str, DeconvolutionResult]:
        """
        Run multiple deconvolution methods.

        Args:
            input_matrix: Input expression matrix (TPM)
            output_dir: Output directory
            methods: List of methods to run (default: all)
            threads: Number of threads
            **kwargs: Method-specific parameters (per method in dict)

        Returns:
            Dictionary mapping method names to results
        """
        if methods is None:
            methods = [
                DeconvolutionMethod.CIBERSORT,
                DeconvolutionMethod.IPS,
                DeconvolutionMethod.ESTIMATE,
                DeconvolutionMethod.MCPCOUNTER,
                DeconvolutionMethod.QUANTISEQ,
                DeconvolutionMethod.EPIC,
            ]

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        results = {}
        for method in methods:
            method_kwargs = kwargs.get(method.value, {})
            output_file = output_dir / f"{method.value}_results.csv"
            result = self.run_deconvolution(
                input_matrix=input_matrix,
                output_file=output_file,
                method=method,
                threads=threads,
                **method_kwargs,
            )
            results[method.value] = result

        return results

    def calculate_signature_scores(
        self,
        input_matrix: Path,
        output_file: Path,
        signature_groups: Union[str, List[str]] = "all",
        method: ScoringMethod = ScoringMethod.INTEGRATION,
        threads: int = 1,
        mini_gene_count: int = 3,
        **kwargs,
    ) -> SignatureScoreResult:
        """
        Calculate signature scores.

        Args:
            input_matrix: Input expression matrix (TPM)
            output_file: Output file path
            signature_groups: Signature group(s) to use
            method: Scoring method
            threads: Number of threads
            mini_gene_count: Minimum genes per signature
            **kwargs: Additional parameters

        Returns:
            SignatureScoreResult
        """
        signature_groups = self.parse_signature_groups(signature_groups)

        cmd = [
            self.iobrpy_cli,
            "calculate_sig_score",
            "--input", str(input_matrix),
            "--output", str(output_file),
            "--signature", ",".join(signature_groups),
            "--method", method.value,
            "--parallel_size", str(threads),
            "--mini_gene_count", str(mini_gene_count),
        ]

        if kwargs.get("adjust_eset"):
            cmd.append("--adjust_eset")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )

            if result.returncode == 0:
                # Try to count signatures
                sig_count = 0
                try:
                    import pandas as pd
                    df = pd.read_csv(output_file)
                    sig_count = len(df.columns) - 1 if 'ID' in df.columns else len(df.columns)
                except Exception:
                    pass

                return SignatureScoreResult(
                    method=method,
                    success=True,
                    output_file=Path(output_file),
                    signature_count=sig_count,
                    message=f"Signature scoring completed",
                )
            else:
                return SignatureScoreResult(
                    method=method,
                    success=False,
                    message=f"Signature scoring failed: {result.stderr}",
                )

        except FileNotFoundError:
            return SignatureScoreResult(
                method=method,
                success=False,
                message=f"iobrpy CLI not found: {self.iobrpy_cli}",
            )
        except Exception as e:
            return SignatureScoreResult(
                method=method,
                success=False,
                message=f"Signature scoring error: {str(e)}",
            )

    def calculate_lr_scores(
        self,
        input_matrix: Path,
        output_file: Path,
        data_type: str = "tpm",
        id_type: str = "ensembl",
        cancer_type: str = "pancan",
        verbose: bool = False,
    ) -> LRResult:
        """
        Calculate ligand-receptor scores.

        Args:
            input_matrix: Input expression matrix (TPM)
            output_file: Output file path
            data_type: Expression type (tpm, fpkm)
            id_type: Gene ID type
            cancer_type: Cancer type
            verbose: Verbose output

        Returns:
            LRResult
        """
        cmd = [
            self.iobrpy_cli,
            "LR_cal",
            "--input", str(input_matrix),
            "--output", str(output_file),
            "--data_type", data_type,
            "--id_type", id_type,
            "--cancer_type", cancer_type,
        ]

        if verbose:
            cmd.append("--verbose")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )

            if result.returncode == 0:
                # Try to count LR pairs
                pair_count = 0
                try:
                    import pandas as pd
                    df = pd.read_csv(output_file)
                    pair_count = len(df)
                except Exception:
                    pass

                return LRResult(
                    success=True,
                    output_file=Path(output_file),
                    pair_count=pair_count,
                    message="LR analysis completed",
                )
            else:
                return LRResult(
                    success=False,
                    message=f"LR analysis failed: {result.stderr}",
                )

        except FileNotFoundError:
            return LRResult(
                success=False,
                message=f"iobrpy CLI not found: {self.iobrpy_cli}",
            )
        except Exception as e:
            return LRResult(
                success=False,
                message=f"LR analysis error: {str(e)}",
            )

    def run_tme_profile(
        self,
        input_matrix: Path,
        output_dir: Path,
        threads: int = 1,
        signature_groups: Union[str, List[str]] = "all",
        deconvolution_methods: Optional[List[DeconvolutionMethod]] = None,
        include_lr: bool = True,
        **kwargs,
    ) -> Dict[str, Any]:
        """
        Run complete TME profile analysis.

        Args:
            input_matrix: Input expression matrix (TPM)
            output_dir: Output directory
            threads: Number of threads
            signature_groups: Signature groups
            deconvolution_methods: Deconvolution methods to run
            include_lr: Include ligand-receptor analysis
            **kwargs: Additional parameters

        Returns:
            Dictionary with all results
        """
        output_dir = Path(output_dir)

        results = {
            'signatures': None,
            'deconvolution': {},
            'lr': None,
        }

        signature_groups = self.parse_signature_groups(signature_groups)
        signatures_dirname = kwargs.get("signatures_dirname", "01-signatures")
        deconvolution_dirname = kwargs.get("deconvolution_dirname", "02-tme")
        lr_dirname = kwargs.get("lr_dirname", "03-LR_cal")

        # 1. Signature scoring
        sig_dir = output_dir / signatures_dirname
        sig_dir.mkdir(parents=True, exist_ok=True)
        sig_file = sig_dir / "signature_scores.csv"

        sig_result = self.calculate_signature_scores(
            input_matrix=input_matrix,
            output_file=sig_file,
            signature_groups=signature_groups,
            threads=threads,
        )
        results['signatures'] = sig_result

        # 2. Deconvolution
        if deconvolution_methods is None:
            deconvolution_methods = [
                DeconvolutionMethod.CIBERSORT,
                DeconvolutionMethod.IPS,
                DeconvolutionMethod.ESTIMATE,
                DeconvolutionMethod.MCPCOUNTER,
                DeconvolutionMethod.QUANTISEQ,
                DeconvolutionMethod.EPIC,
            ]

        deconv_dir = output_dir / deconvolution_dirname
        deconv_dir.mkdir(parents=True, exist_ok=True)

        deconv_results = self.run_all_deconvolution(
            input_matrix=input_matrix,
            output_dir=deconv_dir,
            methods=deconvolution_methods,
            threads=threads,
        )
        results['deconvolution'] = deconv_results

        # 3. Ligand-receptor analysis
        if include_lr:
            lr_dir = output_dir / lr_dirname
            lr_dir.mkdir(parents=True, exist_ok=True)
            lr_file = lr_dir / "lr_scores.csv"

            lr_result = self.calculate_lr_scores(
                input_matrix=input_matrix,
                output_file=lr_file,
                cancer_type=kwargs.get("lr_cancer_type", "pancan"),
            )
            results['lr'] = lr_result

        return results

    def get_available_signatures(self) -> List[str]:
        """
        Get list of available signature groups.

        Returns:
            List of signature group names
        """
        return [
            'go_bp',
            'go_cc',
            'go_mf',
            'signature_collection',
            'signature_tme',
            'signature_sc',
            'signature_tumor',
            'signature_metabolism',
            'kegg',
            'hallmark',
            'reactome',
            'all',
        ]

    def get_deconvolution_cell_types(self, method: DeconvolutionMethod) -> Optional[List[str]]:
        """
        Get cell types for a deconvolution method.

        Args:
            method: Deconvolution method

        Returns:
            List of cell types or None
        """
        cell_types = {
            DeconvolutionMethod.CIBERSORT: [
                'B_cells_naive', 'B_cells_memory', 'Plasma_cells',
                'T_cells_CD8', 'T_cells_CD4_naive', 'T_cells_CD4_memory_resting',
                'T_cells_follicular_helper', 'T_cells_regulatory_Tregs',
                'T_cells_gamma_delta', 'NK_cells_resting', 'NK_cells_activated',
                'Monocytes', 'Macrophages_M0', 'Macrophages_M1',
                'Macrophages_M2', 'Dendritic_cells_resting',
                'Dendritic_cells_activated', 'Mast_cells_resting',
                'Mast_cells_activated', 'Eosinophils', 'Neutrophils',
            ],
            DeconvolutionMethod.QUANTISEQ: [
                'B_cells', 'CD4_T_cells', 'CD8_T_cells', 'Eosinophils',
                'Macrophages', 'Monocytes', 'Neutrophils', 'NK_cells',
                'Tregs', 'Dendritic_cells',
            ],
            DeconvolutionMethod.EPIC: [
                'B_cells', 'CD4_T_cells', 'CD8_T_cells',
                'Endothelial', 'Fibroblasts', 'Macrophages', 'NK_cells',
                'Other_immune', 'Tregs',
            ],
        }
        return cell_types.get(method)

    def parse_deconvolution_result(self, result_file: Path) -> Optional[Dict[str, Dict[str, float]]]:
        """
        Parse a deconvolution result file.

        Args:
            result_file: Path to result file

        Returns:
            Dictionary mapping samples to cell type proportions
        """
        try:
            import pandas as pd
            df = pd.read_csv(result_file)

            if 'ID' not in df.columns:
                return None

            result = {}
            for _, row in df.iterrows():
                sample_id = str(row['ID'])
                result[sample_id] = {
                    col: float(row[col]) for col in df.columns if col != 'ID'
                }
            return result

        except Exception:
            return None
