import copy
from dataclasses import dataclass

import numpy as np
import pandas as pd
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

    def cell_factor_activity(self) -> pd.DataFrame:
        """
        Project each cell onto the protein factors P, after z-scaling it with the same per-protein
        mean and standard deviation (over the pseudobulk group means) as the tensor. The mean over the
        cells of a (treatment, phase) group is then F[t, j, :] with F = G x_1 T x_2 C.
        Returns DataFrame (cells x factors) with treatment and phase columns.
        """
        adata_df = self.adata.to_df()
        pseudobulk_df = adata_df.groupby([self.adata.obs[self.treatment_col], self.adata.obs[self.cc_phase_col]], observed=True).mean()
        scaled = (adata_df - pseudobulk_df.mean()) / pseudobulk_df.std()

        P = self.factors[2]
        df = pd.DataFrame(scaled.values @ P, index=adata_df.index, columns=[f"Factor {f + 1}" for f in range(P.shape[1])])
        df["treatment"] = self.adata.obs[self.treatment_col].astype(str).values
        df["phase"] = self.adata.obs[self.cc_phase_col].astype(str).values
        return df

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

def sign_align_protein_factors(cct_dict):
    """Flip each cell line's protein factors to a common sign across cell lines.

    For factor f, v_f is the leading left singular vector of the (n_proteins x n_cell_lines)
    matrix of loadings [P_1[:, f], ..., P_L[:, f]], oriented to agree with the majority of
    cell lines. Each P_l[:, f] with P_l[:, f] . v_f < 0 is negated together with the core
    slice G_l[:, :, f], which leaves the reconstruction unchanged.

    Returns a dict of aligned (shallow-copied) CCTensor objects and a DataFrame of the signs
    applied (cell lines x factors, -1 = flipped).
    """
    names = list(cct_dict)
    n_factors = cct_dict[names[0]].factors[2].shape[1]
    signs = np.ones((len(names), n_factors))

    for f in range(n_factors):
        V = np.column_stack([cct_dict[name].factors[2][:, f] for name in names])
        dots = np.linalg.svd(V, full_matrices=False)[0][:, 0] @ V
        if dots.sum() < 0:
            dots = -dots
        signs[:, f] = np.where(dots < 0, -1.0, 1.0)

    aligned = {}
    for i, name in enumerate(names):
        cct = copy.copy(cct_dict[name])
        cct.factors = [cct.factors[0], cct.factors[1], cct.factors[2] * signs[i]]
        cct.core = cct.core * signs[i]
        aligned[name] = cct

    columns = [f"Factor {f + 1}" for f in range(n_factors)]
    return aligned, pd.DataFrame(signs.astype(int), index=names, columns=columns)


def cross_cellline_protein_factors(cct_dict):
    """
    Stack protein loadings across cell lines.

    Returns DataFrame of shape (n_protein_factors * n_cell_lines, n_proteins)
    with MultiIndex (cell_line, factor).
    """
    proteins = None
    rows, index_tuples = [], []

    for name, cct in cct_dict.items():
        P = cct.factors[2]
        if proteins is None:
            proteins = list(cct.proteins)
        elif list(cct.proteins) != proteins:
            raise ValueError(f"Protein lists differ: expected {proteins}, got {list(cct.proteins)} for {name}")

        for k in range(P.shape[1]):
            rows.append(P[:, k])
            index_tuples.append((name, f"Factor {k + 1}"))

    index = pd.MultiIndex.from_tuples(index_tuples, names=["cell_line", "factor"])
    return pd.DataFrame(np.vstack(rows), index=index, columns=proteins)


def cross_cellline_treatment_factors(cct_dict):
    """
    Map treatment factors to protein space via core tensor, stack across cell lines.

    For each cell line, computes G x_2 C, averages it over the cell-cycle phases and
    projects to protein space via P. Returns DataFrame of shape
    (n_treatment_factors * n_cell_lines, n_proteins) with MultiIndex (cell_line, factor).
    """
    proteins = None
    rows, index_tuples = [], []

    for name, cct in cct_dict.items():
        P = cct.factors[2]
        if proteins is None:
            proteins = list(cct.proteins)
        elif list(cct.proteins) != proteins:
            raise ValueError(f"Protein lists differ: expected {proteins}, got {list(cct.proteins)} for {name}")

        G_avg = tenalg.mode_dot(cct.core, cct.factors[1], mode=1).mean(axis=1)
        signatures = G_avg @ P.T

        for k in range(signatures.shape[0]):
            rows.append(signatures[k])
            index_tuples.append((name, f"Factor {k + 1}"))

    index = pd.MultiIndex.from_tuples(index_tuples, names=["cell_line", "factor"])
    return pd.DataFrame(np.vstack(rows), index=index, columns=proteins)


def cross_cellline_protein_factors_by_treatment(cct_dict):
    """Scale protein loadings by treatment-specific factor activity, stack across cell lines.

    For each (cell_line, treatment, factor_f), computes F_avg[t, f] * P[:, f] where
    F = G x_1 T x_2 C and F_avg = F.mean(axis=1) averages over cell-cycle phase.

    Returns DataFrame with 3-level MultiIndex (cell_line, treatment, factor) x proteins.
    """
    proteins = None
    rows, index_tuples = [], []

    for name, cct in cct_dict.items():
        P = cct.factors[2]
        if proteins is None:
            proteins = list(cct.proteins)
        elif list(cct.proteins) != proteins:
            raise ValueError(f"Protein lists differ: expected {proteins}, got {list(cct.proteins)} for {name}")

        F = tenalg.multi_mode_dot(cct.core, [cct.factors[0], cct.factors[1]], modes=[0, 1])
        F_avg = F.mean(axis=1)

        for t_idx, treatment in enumerate(cct.treatments):
            for f in range(P.shape[1]):
                rows.append(F_avg[t_idx, f] * P[:, f])
                index_tuples.append((name, treatment, f"Factor {f + 1}"))

    index = pd.MultiIndex.from_tuples(index_tuples, names=["cell_line", "treatment", "factor"])
    return pd.DataFrame(np.vstack(rows), index=index, columns=proteins)


def cross_cellline_phase_profiles(cct_dict):
    """Phase-resolved factor activations across cell lines.

    For each (cell_line, treatment, factor_f), stores F[t, :, f] where
    F = G x_1 T x_2 C. Each row is a 3-dimensional vector giving the
    factor's activation in each cell-cycle phase (not averaged).

    Returns DataFrame with 3-level MultiIndex (cell_line, treatment, factor) x phases.
    """
    phases = None
    rows, index_tuples = [], []

    for name, cct in cct_dict.items():
        if phases is None:
            phases = list(cct.cc_phases)

        F = tenalg.multi_mode_dot(cct.core, [cct.factors[0], cct.factors[1]], modes=[0, 1])

        for t_idx, treatment in enumerate(cct.treatments):
            for f in range(F.shape[2]):
                rows.append(F[t_idx, :, f])
                index_tuples.append((name, treatment, f"Factor {f + 1}"))

    index = pd.MultiIndex.from_tuples(index_tuples, names=["cell_line", "treatment", "factor"])
    return pd.DataFrame(np.vstack(rows), index=index, columns=phases)


def cross_cellline_cell_activity(cct_dict, random_state=0):
    """Per-cell factor activity (CCTensor.cell_factor_activity) stacked across cell lines.

    For each treatment, every cell line contributes the same number of randomly drawn cells: the
    count of the smallest line under that treatment. Factors are assumed sign-aligned
    (sign_align_protein_factors).

    Returns DataFrame (cells x factors) with treatment, phase and cell_line columns.
    """
    df = pd.concat(
        [cct.cell_factor_activity().assign(cell_line=name) for name, cct in cct_dict.items()],
        ignore_index=True,
    )
    n = df.groupby(["treatment", "cell_line"]).size().groupby(level="treatment").min()
    df = df.sample(frac=1, random_state=random_state)
    keep = df.groupby(["treatment", "cell_line"]).cumcount() < df["treatment"].map(n)
    return df[keep].sort_index()


def cross_cellline_protein_influence(cct_dict):
    """Phase-resolved protein influence of each latent protein-factor, stacked across cell lines.

    For each (cell_line, treatment, phase, factor_f), computes F[t, ph, f] * P[:, f] where
    F = G x_1 T x_2 C. Averaging over phase recovers cross_cellline_protein_factors_by_treatment.
    The product is invariant to the sign ambiguity of the HOSVD loadings, so rows can be pooled
    across cell lines.

    Returns DataFrame with 4-level MultiIndex (cell_line, treatment, phase, factor) x proteins.
    """
    proteins = None
    rows, index_tuples = [], []

    for name, cct in cct_dict.items():
        P = cct.factors[2]
        if proteins is None:
            proteins = list(cct.proteins)
        elif list(cct.proteins) != proteins:
            raise ValueError(f"Protein lists differ: expected {proteins}, got {list(cct.proteins)} for {name}")

        F = tenalg.multi_mode_dot(cct.core, [cct.factors[0], cct.factors[1]], modes=[0, 1])

        for t_idx, treatment in enumerate(cct.treatments):
            for ph_idx, phase in enumerate(cct.cc_phases):
                for f in range(P.shape[1]):
                    rows.append(F[t_idx, ph_idx, f] * P[:, f])
                    index_tuples.append((name, treatment, phase, f"Factor {f + 1}"))

    index = pd.MultiIndex.from_tuples(index_tuples, names=["cell_line", "treatment", "phase", "factor"])
    return pd.DataFrame(np.vstack(rows), index=index, columns=proteins)
