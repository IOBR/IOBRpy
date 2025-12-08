import numpy as np
import pandas as pd
import warnings

from .gibbs import GibbsSampler
from .joint_post import JointPost
from . import process_input
from . import references as rf
from . import optim

class Prism:
    def __init__(self, phi_cellState, phi_cellType, map_, key, mixture):
        self.phi_cellState = phi_cellState
        self.phi_cellType = phi_cellType
        self.map = map_
        self.key = key
        self.mixture = mixture
    
    def valid_opt_control(control):
        ctrl = {'maxit': 100000, 
                'maximize': False, 
                'trace': 0, 
                'eps': 1e-07,
                'dowarn': True, 
                'tol': 0, 
                'maxNA': 500, 
                'n.cores': 1, 
                'optimizer': 'MAP', 
                'sigma': 2}
        
        namc = list(control.keys())

        for name in namc:
            if name in ctrl.keys():
                ctrl[name] = control[name]
            else:
                raise ValueError("Unknown names in opt.control: {}".format(name))

        if ctrl['optimizer'] not in ["MAP", "MLE"]:
            raise ValueError("unknown names of optimizer: " + ctrl['optimizer'])
        
        if ctrl['optimizer'] == "MAP":
            if not isinstance(ctrl['sigma'], (int, float)):
                raise ValueError("sigma needs to be a numeric variable")
            else:
                if ctrl['sigma'] < 0:
                    raise ValueError("sigma needs to be positive")
        return ctrl


    def valid_gibbs_control(control):
        ctrl = {'chain.length': 1000,
                'burn.in': 500,
                'thinning': 2,
                'n.cores': 1,
                'seed': 123,
                'alpha': 1,
                'fast.multinomial': False,
                'rng.backend': 'generator'}

        namc = list(control.keys())
       
        for name in namc:
            if name in ctrl.keys():
                ctrl[name] = control[name]
            else:
                raise ValueError("Unknown names in opt.control: {}".format(name))

        if ctrl['alpha'] < 0:
            raise ValueError("alpha needs to be positive")

        if not isinstance(ctrl['fast.multinomial'], bool):
            raise ValueError("fast.multinomial must be a boolean flag")

        if ctrl['rng.backend'].lower() not in ['generator', 'randomstate']:
            raise ValueError("rng.backend must be either 'generator' or 'randomstate'")

        ctrl['rng.backend'] = ctrl['rng.backend'].lower()

        return ctrl


    def new(reference, input_type, cell_type_labels, cell_state_labels, \
            key, mixture, outlier_cut=0.01, outlier_fraction=0.1, pseudo_min=1E-8):
        
        if cell_state_labels is None:
            cell_state_labels = cell_type_labels

        print("number of cells in each cell state")
        print(pd.Series(cell_state_labels).value_counts().sort_values(ascending = False))
        if np.min(pd.Series(cell_state_labels).value_counts()) < 20:
            print("recommend to have sufficient number of cells in each cell state")
        
        if key is None:
            print("No tumor reference is speficied. Reference cell types are treated equally.")
        
        if len(cell_type_labels) != len(cell_state_labels):
            raise ValueError("Error: length of cell.type.labels and cell.state.labels do not match!")
        if len(cell_type_labels) != reference.shape[0]:
            raise ValueError("Error: length of cell.type.labels and nrow(reference) do not match!")
        
        type_to_state_mat = pd.DataFrame({"cell.type.labels": cell_type_labels, "cell.state.labels": cell_state_labels})
        type_to_state_mat = type_to_state_mat.drop_duplicates()
        if type_to_state_mat["cell.state.labels"].value_counts().max() > 1:
            raise ValueError("Error: one or more cell states belong to multiple cell types!")
        if len(pd.unique(cell_type_labels)) > len(pd.unique(cell_state_labels)):
            raise ValueError("Error: more cell types than states!")
        
        if not isinstance(mixture, pd.DataFrame):
            mixture = pd.DataFrame(mixture)
        if not isinstance(reference, pd.DataFrame):
            reference = pd.DataFrame(reference)

        if mixture.shape[1] == 1:
            mixture = mixture.T
            mixture = mixture.rename(index = {0: ["mixture-1"]})

        if mixture.index.equals(pd.RangeIndex(mixture.shape[0])):
            mixture.index = [f"mixture-{i}" for i in range(len(mixture))]
        
        process_input.validate_input(reference)
        process_input.validate_input(mixture)

        input_type = input_type.lower()
        if input_type == "count.matrix":
            # Keep genes expressed in at least one cell to mirror R BayesPrism's
            # count-matrix filtering behaviour.
            expressed_mask = np.sum(reference > 0, axis=0) > 0
            if expressed_mask.sum() < reference.shape[1]:
                print(f"Removing {reference.shape[1] - expressed_mask.sum()} genes not expressed in any cell")
            reference = reference.loc[:, expressed_mask]
            pseudo_min_to_use = pseudo_min
        elif input_type in ["tpm", "gep"]:
            # Collapsed/normalized profiles (e.g., TPM) only support removing
            # genes with zero abundance across all cell states/types and do not
            # add pseudo count during normalization.
            pseudo_min_to_use = 0
            if pseudo_min != 0:
                warnings.warn(
                    "pseudo_min is ignored for normalized input types (TPM/GEP); using 0 instead."
                )
            reference = reference.loc[:, np.sum(reference, axis=0) > 0]
        else:
            raise ValueError(
                "Unsupported input_type. Please specify one of ['count.matrix', 'TPM', 'GEP']."
            )

        mixture = process_input.filter_bulk_outlier(mixture, outlier_cut, outlier_fraction)

        _, gene_index, _ = np.intersect1d(reference.columns, mixture.columns, return_indices = True)
        gene_index.sort()
        gene_shared = [list(reference.columns)[i] for i in gene_index]
        if len(gene_shared) == 0:
            raise ValueError("Error: gene names of reference and mixture do not match!")
        if len(gene_shared) < 100:
            print("Warning: very few gene from reference and mixture match! Please double check your gene names.")
        
        ref_cs = process_input.collapse(reference, cell_state_labels)
        ref_ct = process_input.collapse(reference, cell_type_labels)
        
        print("Aligning reference and mixture...")
        ref_ct = ref_ct.loc[:, gene_shared]
        ref_cs = ref_cs.loc[:, gene_shared]
        mixture = mixture.loc[:, gene_shared]
        
        print("Normalizing reference...")
        ref_cs = process_input.norm_to_one(ref_cs, pseudo_min_to_use)
        ref_ct = process_input.norm_to_one(ref_ct, pseudo_min_to_use)
        
        map_ = {cell_type: list(set(cell_state_labels[i] for i, ct in enumerate(cell_type_labels) if ct == cell_type)) for cell_type in ref_ct.index}

        return Prism(
            rf.RefPhi(ref_cs, pseudo_min_to_use),
            rf.RefPhi(ref_ct, pseudo_min_to_use),
            map_,
            key,
            mixture
        )

    def run(self, n_cores = 1, update_gibbs = True, gibbs_control = {}, opt_control = {}):
        
        if 'n.cores' not in gibbs_control:
            gibbs_control['n.cores'] = n_cores
        if 'n.cores' not in opt_control:
            opt_control['n.cores'] = n_cores
        
        assert isinstance(update_gibbs, bool)
        assert isinstance(n_cores, int)
        
        opt_control = Prism.valid_opt_control(opt_control)
        gibbs_control = Prism.valid_gibbs_control(gibbs_control)
        
        if self.phi_cellState.pseudo_min == 0:
            gibbs_control['alpha'] = max(1, gibbs_control['alpha'])
        
        gibbsSampler_ini_cs = GibbsSampler(reference = self.phi_cellState, 
                                           X = self.mixture, 
                                           gibbs_control = gibbs_control)

        jointPost_ini_cs = gibbsSampler_ini_cs.run(final = False)
        print("Now Merging...")
        jointPost_ini_ct = jointPost_ini_cs.merge_K(map_ = self.map)

        if not update_gibbs:
            bp = BayesPrism(prism = self, 
                                posterior_initial_cellState = jointPost_ini_cs,
                                posterior_initial_cellType = jointPost_ini_ct,
                                control_param = {'gibbs.control': gibbs_control,
                                                 'opt.control': opt_control,
                                                 'update.gibbs': update_gibbs})
            return bp
        else:
            psi = optim.update_reference(Z = jointPost_ini_ct.Z, 
                                         phi_prime = self.phi_cellType, 
                                         map = self.map,
                                         key = self.key, 
                                         opt_control = opt_control)
            
            gibbsSampler_update = GibbsSampler(reference = psi, 
                                               X = self.mixture, 
                                               gibbs_control = gibbs_control)
            
            theta_f = gibbsSampler_update.run(final = True)

            bp = BayesPrism(prism = self, 
                                posterior_initial_cellState = jointPost_ini_cs,
                                posterior_initial_cellType = jointPost_ini_ct, 
                                control_param={'gibbs.control': gibbs_control,
                                               'opt.control': opt_control,
                                               'update.gibbs': update_gibbs},
                                reference_update = psi,
                                posterior_theta_f = theta_f)
            return bp




class BayesPrism:
    def __init__(self, prism, posterior_initial_cellState, posterior_initial_cellType,
                 control_param, reference_update = None, posterior_theta_f = None):
        self.prism = prism
        self.posterior_initial_cellState = posterior_initial_cellState
        self.posterior_initial_cellType = posterior_initial_cellType
        self.control_param = control_param
        self.reference_update = reference_update
        self.posterior_theta_f = posterior_theta_f


    def update_theta(self, gibbs_control = {}, opt_control = {}):

        gibbs_control_bp = self.control_param['gibbs.control']
        opt_control_bp = self.control_param['opt.control']
        
        if gibbs_control:
            for key, value in gibbs_control.items():
                gibbs_control_bp[key] = value
    
        if opt_control:
            for key, value in opt_control.items():
                opt_control_bp[key] = value

        opt_control = Prism.valid_opt_control(opt_control_bp)
        gibbs_control = Prism.valid_gibbs_control(gibbs_control_bp)
        
        psi = optim.update_reference(Z = self.posterior_initial_cellType.Z, 
                                     phi_prime = self.prism.phi_cellType,
                                     map = self.prism.map, 
                                     key = self.prism.key, 
                                     opt_control = opt_control)
        
        gibbsSampler_update = GibbsSampler(reference = psi, 
                                           X = self.prism.mixture, 
                                           gibbs_control = gibbs_control)
        
        theta_f = gibbsSampler_update.run(final = True)
        
        bp_updated = BayesPrism(prism = self.prism, 
                                posterior_initial_cellState = self.posterior_initial_cellState,
                                posterior_initial_cellType = self.posterior_initial_cellType, 
                                control_param={'gibbs.control': gibbs_control, 
                                               'opt.control': opt_control, 
                                               'update.gibbs': True},
                                reference_update = psi,
                                posterior_theta_f = theta_f)

        return bp_updated