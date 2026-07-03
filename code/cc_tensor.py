from dataclasses import dataclass

import numpy as np
import tensorly as tl
from tensorly import tenalg

import matplotlib.pyplot as plt


@dataclass
class CCTensor:
    """Class for manipulating cell-cycle tensors."""
    adata: any = None
    
    treatment_col: str = "drug_name"
    cc_phase_col: str = "GMM_phase_labels"
    
    treatments: list = None
    cc_phases: list = None
    proteins: list = None
    
    tensor: any = None
    core: any = None
    factors: list = None


    def __init__(
            self, 
            adata: any, 
            treatment_col: str = "drug_name", cc_phase_col: str = "GMM_phase_labels", 
            treatments: list = None, cc_phases: list = None, proteins: list = None
        ) -> None:
        self.adata = adata
        
        self.treatment_col = treatment_col
        self.cc_phase_col = cc_phase_col
        
        self.treatments = treatments
        self.cc_phases = cc_phases
        self.proteins = proteins

    def pseudobulk(self) -> any:
        """
        Create a pseudobulk tensor by aggregating cells by treatment and cell cycle phase.
        Returns z-scaled tensor of shape (n_treatments, n_phases, n_proteins).
        """
        # Convert AnnData to DataFrame 
        adata_df = self.adata.to_df()
        
        # Group by treatment_col and cc_phase_col, then get mean of each group
        pseudobulk_df = adata_df.groupby([self.adata.obs[self.treatment_col], self.adata.obs[self.cc_phase_col]], observed=True).mean()
        
        # z-scale 
        pseudobulk_df_scales = (pseudobulk_df - pseudobulk_df.mean()) / pseudobulk_df.std()
        
        # Reshape to (n_treatments, n_phases, n_proteins)
        return_tensor = pseudobulk_df_scales.values.reshape(
            len(self.adata.obs[self.treatment_col].cat.categories), 
            len(self.adata.obs[self.cc_phase_col].cat.categories), 
            -1
        )
        
        self.tensor = return_tensor
        return return_tensor
    
    def HOSVD(self, truncation_rank=None) -> tuple:
        """
        Compute the Higher-Order Singular Value Decomposition (HOSVD) of the tensor.
        Returns the core tensor and factor matrices.
        """
        if self.tensor is None:
            raise ValueError("Tensor is not initialized. Please run pseudobulk() first.")
        
        dims = self.tensor.shape

        if truncation_rank is None:
            truncation_rank = dims
        
        factors = []
        for i, _ in enumerate(dims): 
            matrix = tl.unfold(self.tensor, mode=i)
            
            rank = np.linalg.matrix_rank(matrix)
            factor = np.linalg.svd(matrix, full_matrices=False)[0][:, :rank]
            
            factors.append(factor[:, :truncation_rank[i]])
        
        factors_transposed = [f.T for f in factors]
        
        core = tenalg.multi_mode_dot(self.tensor, factors_transposed)
            
        return core, factors
    
    def set_tensor_decompositions(self, truncation_rank=None) -> None:
        """
        Compute and store the HOSVD of the tensor.
        """
        core, factors = self.HOSVD(truncation_rank=truncation_rank)
        self.core = core
        self.factors = factors
    
    def determine_rank(self) -> int:
        # Iterate over protein-ranks and record reconstruction errors
        err_dict = {}
        for i in range(len(self.proteins)):
            truncation_rank = [self.tensor.shape[0], self.tensor.shape[1], i+1]
            core, factors = self.HOSVD(truncation_rank=truncation_rank)
            
            recon = tenalg.multi_mode_dot(core, factors)
            err = np.linalg.norm(recon - self.tensor) / np.linalg.norm(self.tensor)
            
            err_dict[i+1] = err
            
            
        # We take as an approximation to the optimal rank the point where the difference between the 
        # reconstruction error and a linear interpolation between the first and last points is maximized.
        start_err = err_dict[1]
        end_err = err_dict[len(err_dict)]
        steps = len(err_dict)

        interpol_line = {i: (start_err * (((steps-1) - (i-1))/(steps-1)) + end_err * ((i-1)/(steps-1))) for i in range(1, steps+1)}
        diff = {i: interpol_line[i] - err_dict[i] for i in range(1, steps+1)}

        n_protein_factors = max(diff, key=diff.get)
        print(f"Maximum difference at step {n_protein_factors}: {diff[n_protein_factors]}")
        
        
        # Plotting the reconstruction error and the interpolated line
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(1, 1, 1)

        ax.set_ylabel("Reconstruction Error")
        ax.set_xlabel("Protein Rank")
        ax.set_xticks([k for k in err_dict if k % 4 == 0])

        ax.spines[["top", "right"]].set_visible(False)

        ax.plot(list(err_dict.keys()), list(err_dict.values()), marker='o')
        ax.plot(list(interpol_line.keys()), list(interpol_line.values()), alpha=0.7, marker='x', linestyle='--')
        ax.plot(list(diff.keys()), list(diff.values()), alpha=0.7, marker='^', linestyle='--')
        ax.vlines(n_protein_factors, 0, err_dict[1], colors='red', alpha=0.25, linestyle='--')

        plt.legend(["Reconstruction Error", "Interpolated", "Difference"])
        plt.title("Reconstruction Error vs. Protein Rank")
        plt.show()
        
        
        # Set tensor decomposition and return the number of protein factors
        self.set_tensor_decompositions(truncation_rank=[self.tensor.shape[0], self.tensor.shape[1], n_protein_factors])
        return n_protein_factors
        
        