# ./workflow/__init__.py

from .prepare_salmon import prepare_salmon_tpm as prepare_salmon_tpm_main
from .count2tpm import count2tpm as count2tpm_main
from .anno_eset import main as anno_eset_main
from .calculate_sig_score import calculate_sig_score as calculate_sig_score_main
from .cibersort import cibersort as cibersort_main
from .IPS import main as IPS_main
from .estimate import estimate_score as estimate_score_main
from .mcpcounter import MCPcounter_estimate as MCPcounter_estimate_main
from .mcpcounter import preprocess_input as preprocess_input_main
from .quantiseq import main as quantiseq_main
from .epic import main as epic_main
from .deside import main as deside_main
from .tme_cluster import main as tme_cluster_main
from .LR_cal import main as LR_cal_main
from .nmf import main as nmf_main