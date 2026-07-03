# %% Imports and constants
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import scanpy as sc
import phate
import tensorly as tl
from tensorly import tenalg
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import gaussian_kde
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

POS = "#b2182b"
NEG = "#2166ac"
CMAP = "RdBu_r"
PHASE_COLORS = {"G1": "#2ca02c", "S": "#ff7f0e", "G2": "#d62728"}

_RC = {
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "figure.facecolor": "white",
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "image.cmap": CMAP,
}


# %% CCTensorViz
class CCTensorViz:

    def __init__(self, cct, n_sub=10_000, n_per=3000, dot_size=5, random_state=0):
        if cct.core is None or cct.factors is None:
            raise ValueError("CCTensor must be fitted (call set_tensor_decompositions or determine_rank first).")
        if cct.tensor is None:
            raise ValueError("CCTensor must have a tensor (call pseudobulk first).")

        self.cct = cct
        self.n_sub = n_sub
        self.n_per = n_per
        self.dot_size = dot_size
        self.random_state = random_state

        self.treatments = np.asarray(cct.treatments)
        self.cc_phases = np.asarray(cct.cc_phases)
        self.proteins = np.asarray(cct.proteins)

        self.treatment_loadings = cct.factors[0]
        self.phase_loadings = cct.factors[1]
        self.protein_loadings = cct.factors[2]

        self.n_factors = self.protein_loadings.shape[1]
        self.factor_names = [f"Factor {i + 1}" for i in range(self.n_factors)]

        self.F = tenalg.multi_mode_dot(
            cct.core,
            [self.treatment_loadings, self.phase_loadings],
            modes=[0, 1],
        )

        X_full = np.asarray(cct.adata.X)
        self._mu = X_full.mean(0)
        self._sd = np.where(X_full.std(0) == 0, 1.0, X_full.std(0))

        self._global_phate = None
        self._adata_sub = None
        self._global_scores = None
        self._per_treatment = None
        self._per_treatment_dpt_done = False

        plt.rcParams.update(_RC)

    # ------------------------------------------------------------------ #
    # Private helpers
    # ------------------------------------------------------------------ #

    @staticmethod
    def _grid_shape(n, max_cols=4):
        cols = min(n, max_cols)
        rows = math.ceil(n / cols)
        return rows, cols

    def _annotate_heatmap(self, ax, matrix, vmax, text_size=8):
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                v = matrix[i, j]
                ax.text(
                    j, i, f"{v:.2f}", ha="center", va="center", fontsize=text_size,
                    color="white" if abs(v) > 0.6 * vmax else "black",
                )

    def _embed_factors(self, X):
        Xz = (np.asarray(X) - self._mu) / self._sd
        return Xz, Xz @ self.protein_loadings

    def _ensure_global_phate(self):
        if self._global_phate is not None:
            return
        print("Computing global PHATE embedding...")
        adata = self.cct.adata
        self._adata_sub = sc.pp.subsample(
            adata, n_obs=min(self.n_sub, adata.n_obs),
            copy=True, random_state=self.random_state,
        )
        phate_op = phate.PHATE(random_state=self.random_state)
        self._global_phate = phate_op.fit_transform(self._adata_sub.X)
        _, self._global_scores = self._embed_factors(self._adata_sub.X)

    def _ensure_per_treatment(self):
        if self._per_treatment is not None:
            return
        print("Computing per-treatment PHATE embeddings...")
        adata = self.cct.adata
        tcol = self.cct.treatment_col
        pcol = self.cct.cc_phase_col

        self._per_treatment = {}
        for t in self.treatments:
            mask = adata.obs[tcol] == t
            sub = sc.pp.subsample(
                adata[mask],
                n_obs=min(self.n_per, int(mask.sum())),
                copy=True, random_state=self.random_state,
            )
            Xz, embed = self._embed_factors(sub.X)
            norms = np.linalg.norm(Xz, axis=1, keepdims=True)

            self._per_treatment[t] = {
                "adata": sub,
                "Xz": Xz,
                "embed": embed,
                "phate_coords": phate.PHATE(random_state=self.random_state, verbose=0).fit_transform(Xz),
                "cos_dist": 1 - embed / norms,
                "phase": sub.obs[pcol].astype(str).values,
            }

    def _ensure_per_treatment_dpt(self):
        if self._per_treatment_dpt_done:
            return
        self._ensure_per_treatment()
        print("Computing per-treatment diffusion pseudotime...")
        pcol = self.cct.cc_phase_col

        for t in self.treatments:
            pt = self._per_treatment[t]
            sub = pt["adata"].copy()
            sub.X = pt["Xz"].copy()

            sc.pp.neighbors(sub, random_state=self.random_state)  # Build KNN graph
            sc.tl.diffmap(sub, n_comps=10)  # Transition matrix of random walks on the KNN graph

            g1_mask = sub.obs[pcol].astype(str).values == "G1"
            if g1_mask.any():
                dc1 = sub.obsm["X_diffmap"][:, 0]
                g1_indices = np.flatnonzero(g1_mask)
                sub.uns["iroot"] = int(g1_indices[np.argmin(dc1[g1_mask])])
            else:
                sub.uns["iroot"] = 0

            sc.tl.dpt(sub, n_dcs=10)
            pt["dpt_pseudotime"] = sub.obs["dpt_pseudotime"].values

        self._per_treatment_dpt_done = True

    @staticmethod
    def _add_phase_contours(ax, coords, phases, enclose_frac=0.80,
                            min_cells=15, linewidth=1.5, pad_frac=0.1):
        xmin, xmax = coords[:, 0].min(), coords[:, 0].max()
        ymin, ymax = coords[:, 1].min(), coords[:, 1].max()
        dx = (xmax - xmin) * pad_frac
        dy = (ymax - ymin) * pad_frac
        xx, yy = np.mgrid[
            xmin - dx : xmax + dx : 100j,
            ymin - dy : ymax + dy : 100j,
        ]
        grid_pts = np.vstack([xx.ravel(), yy.ravel()])

        for phase, color in PHASE_COLORS.items():
            mask = phases == phase
            if mask.sum() < min_cells:
                continue
            pts = coords[mask].T
            try:
                kde = gaussian_kde(pts)
            except np.linalg.LinAlgError:
                continue
            densities_at_data = kde(pts)
            level = np.percentile(densities_at_data, 100 * (1 - enclose_frac))
            if level <= 0:
                continue
            Z = kde(grid_pts).reshape(xx.shape)
            ax.contour(xx, yy, Z, levels=[level], colors=[color],
                       linewidths=linewidth, alpha=0.9)

    # ------------------------------------------------------------------ #
    # Tensor-only visualizations (no PHATE)
    # ------------------------------------------------------------------ #

    def protein_loadings_heatmap(self, annotate=True, text_size=8):
        fig, ax = plt.subplots(figsize=(max(4, self.n_factors * 2), max(6, len(self.proteins) * 0.55)))

        vmax = np.abs(self.protein_loadings).max() or 1.0
        im = ax.imshow(self.protein_loadings, cmap=CMAP, vmin=-vmax, vmax=vmax, aspect="auto")

        ax.set_xticks(np.arange(self.n_factors))
        ax.set_xticklabels([f"{i+1}" for i in range(self.n_factors)])
        ax.set_xlabel("Latent Protein-Factors", fontsize=10)

        ax.set_yticks(np.arange(len(self.proteins)))
        ax.set_yticklabels(self.proteins, fontsize=8)
        ax.set_ylabel("Protein", fontsize=10)

        if annotate:
            self._annotate_heatmap(ax, self.protein_loadings, vmax, text_size)

        fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.4, pad=0.05)
        ax.set_title("Protein Loadings", fontsize=13)
        fig.tight_layout()
        return fig, ax

    def treatment_loadings_heatmap(self, annotate=True, text_size=8):
        T = self.treatment_loadings
        n_t_factors = T.shape[1]

        fig, ax = plt.subplots(figsize=(max(4, n_t_factors * 1.5), max(6, len(self.treatments) * 0.55)))

        vmax = np.abs(T).max() or 1.0
        im = ax.imshow(T, cmap=CMAP, vmin=-vmax, vmax=vmax, aspect="auto")

        ax.set_xticks(np.arange(n_t_factors))
        ax.set_xticklabels([f"{i+1}" for i in range(n_t_factors)])
        ax.set_xlabel("Latent Treatment-Factors", fontsize=10)

        ax.set_yticks(np.arange(len(self.treatments)))
        ax.set_yticklabels(self.treatments, fontsize=8)
        ax.set_ylabel("Treatment", fontsize=10)

        if annotate:
            self._annotate_heatmap(ax, T, vmax, text_size)

        fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.4, pad=0.05)
        ax.set_title("Treatment Loadings", fontsize=13)
        fig.tight_layout()
        return fig, ax

    def protein_loadings_lollipop(self):
        n = self.n_factors
        vmax = np.abs(self.protein_loadings).max() or 1.0

        rows, cols = self._grid_shape(n)
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, max(6, len(self.proteins) * 0.3)))
        axes = np.atleast_1d(axes).ravel()

        for k, ax in enumerate(axes):
            if k >= n:
                ax.axis("off")
                continue

            order = np.argsort(self.protein_loadings[:, k])
            vals = self.protein_loadings[order, k]
            y = np.arange(len(vals))
            colors = np.where(vals >= 0, POS, NEG)

            ax.hlines(y, 0, vals, color=colors, lw=1.5)
            ax.scatter(vals, y, color=colors, s=30, zorder=3)

            ax.set_yticks(y)
            ax.set_yticklabels(self.proteins[order], fontsize=7)
            ax.set_ylim(-0.5, len(vals) - 0.5)
            ax.set_xlim(-1.05 * vmax, 1.05 * vmax)
            ax.axvline(0, color="0.7", lw=0.8, zorder=0)
            ax.set_title(f"Factor {k + 1}", fontsize=10)

        sns.despine(fig=fig)
        fig.suptitle("Latent protein-factor composition", y=1.0, fontsize=13)
        fig.tight_layout()
        return fig, axes

    def gram_clustermap(self, mode="protein", method="average", annotate=True):
        if mode == "protein":
            M = self.protein_loadings
            labels = list(self.proteins)
            title = "Protein cosine similarity (latent factor space)"
        else:
            M = self.treatment_loadings
            labels = list(self.treatments)
            title = "Treatment cosine similarity (latent factor space)"

        norms = np.linalg.norm(M, axis=1, keepdims=True)
        norms = np.where(norms == 0, 1.0, norms)
        M_norm = M / norms

        gram = M_norm @ M_norm.T
        dist = (1 - gram) / 2
        np.fill_diagonal(dist, 0)

        condensed = squareform(dist, checks=False)
        Z = linkage(condensed, method=method)

        g = sns.clustermap(
            gram,
            row_linkage=Z,
            col_linkage=Z,
            cmap=CMAP,
            center=0,
            vmin=-1,
            vmax=1,
            annot=annotate,
            fmt=".2f",
            annot_kws={"size": 7},
            xticklabels=labels,
            yticklabels=labels,
            figsize=(max(6, len(labels) * 0.6), max(6, len(labels) * 0.6)),
            dendrogram_ratio=0.15,
            cbar_kws={"label": "Cosine similarity"},
        )

        g.cax.set_visible(False)
        g.ax_col_dendrogram.set_title(title, fontsize=12, pad=10)
        g.ax_heatmap.tick_params(axis="both", labelsize=8)
        plt.tight_layout()

        return g

    def F_slice_montage(self, annotate=True):
        n = self.n_factors
        cols = min(n, 2)
        rows = math.ceil(n / cols)
        vmax = np.abs(self.F).max() or 1.0

        fig, axes = plt.subplots(rows, cols, figsize=(3 * n, 4 * rows), squeeze=False)
        flat = axes.ravel()

        for f in range(n):
            ax = flat[f]
            im = ax.imshow(self.F[:, :, f], cmap=CMAP, vmin=-vmax, vmax=vmax, aspect="auto")

            ax.set_title(f"Latent Protein-Factor {f + 1}", fontsize=10)
            ax.set_xticks(np.arange(len(self.cc_phases)))
            ax.set_xticklabels(self.cc_phases, fontsize=8)
            ax.set_xlabel("Phase", fontsize=8)
            ax.set_yticks(np.arange(len(self.treatments)))
            show_yticks = (f % cols == 0)
            ax.set_yticklabels(self.treatments if show_yticks else [], fontsize=8)

            if annotate:
                self._annotate_heatmap(ax, self.F[:, :, f], vmax, text_size=6)

        for idx in range(n, len(flat)):
            flat[idx].axis("off")

        fig.suptitle("F -- Treatment x Latent Protein-Factors by Phase", fontsize=13, y=0.94, x=0.42)
        fig.tight_layout(rect=[0, 0, 0.88, 0.96])
        fig.colorbar(im, ax=flat[:n].tolist(), orientation="vertical", shrink=0.8, pad=0.04)
        return fig, axes

    def F_slice_montage_normalized(self, control="Control", annotate=True):
        matches = np.flatnonzero(self.treatments == control)
        if len(matches) == 0:
            raise ValueError(f"Control treatment '{control}' not found in {list(self.treatments)}.")
        ctrl_idx = matches[0]

        F_delta = self.F - self.F[ctrl_idx : ctrl_idx + 1, :, :]

        n = self.n_factors
        cols = min(n, 2)
        rows = math.ceil(n / cols)
        vmax = np.abs(F_delta).max() or 1.0

        fig, axes = plt.subplots(rows, cols, figsize=(3 * n, 4 * rows), squeeze=False)
        flat = axes.ravel()

        for f in range(n):
            ax = flat[f]
            im = ax.imshow(F_delta[:, :, f], cmap=CMAP, vmin=-vmax, vmax=vmax, aspect="auto")

            ax.set_title(f"Latent Protein-Factor {f + 1}", fontsize=10)
            ax.set_xticks(np.arange(len(self.cc_phases)))
            ax.set_xticklabels(self.cc_phases, fontsize=8)
            ax.set_xlabel("Phase", fontsize=8)
            ax.set_yticks(np.arange(len(self.treatments)))
            show_yticks = (f % cols == 0)
            ax.set_yticklabels(self.treatments if show_yticks else [], fontsize=8)

            if annotate:
                self._annotate_heatmap(ax, F_delta[:, :, f], vmax, text_size=6)

        for idx in range(n, len(flat)):
            flat[idx].axis("off")

        fig.suptitle(
            f"F (Control-Normalized) -- Treatment Effect on Latent Protein-Factors by Phase",
            fontsize=13, y=0.94, x=0.42,
        )
        fig.tight_layout(rect=[0, 0, 0.88, 0.96])
        fig.colorbar(im, ax=flat[:n].tolist(), orientation="vertical", shrink=0.8, pad=0.04)
        return fig, axes

    def singular_value_spectrum(self):
        mode_labels = ["Treatment", "Cell-cycle phase", "Protein"]
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))

        for i, (ax, label) in enumerate(zip(axes, mode_labels)):
            matrix = tl.unfold(self.cct.tensor, mode=i)
            sv = np.linalg.svd(matrix, compute_uv=False)

            total = (sv ** 2).sum()
            cumvar = np.cumsum(sv ** 2) / total

            ax.stem(np.arange(1, len(sv) + 1), sv, linefmt="C0-", markerfmt="C0o", basefmt="C0-")
            ax.set_xlabel("Component")
            ax.set_ylabel("Singular value")
            ax.set_title(label)
            ax.set_xticks(np.arange(1, len(sv) + 1))

            ax2 = ax.twinx()
            ax2.plot(np.arange(1, len(sv) + 1), cumvar, "C1x--", alpha=0.7)
            ax2.set_ylabel("Cumulative variance", color="C1")
            ax2.set_ylim(0, 1.05)
            ax2.tick_params(axis="y", labelcolor="C1")

        fig.suptitle("Singular value spectrum per mode", fontsize=13, y=1.02)
        fig.tight_layout()
        return fig, axes

    def reconstruction_error(self, max_rank=None):
        if max_rank is None:
            max_rank = len(self.proteins)

        tensor = self.cct.tensor
        err_dict = {}
        for k in range(1, max_rank + 1):
            tr = [tensor.shape[0], tensor.shape[1], k]
            core, factors = self.cct.HOSVD(truncation_rank=tr)
            recon = tenalg.multi_mode_dot(core, factors)
            err_dict[k] = np.linalg.norm(recon - tensor) / np.linalg.norm(tensor)

        ranks = list(err_dict.keys())
        errs = list(err_dict.values())

        start, end = errs[0], errs[-1]
        n = len(ranks)
        interpol = [start + (end - start) * i / (n - 1) for i in range(n)]
        diff = [interpol[i] - errs[i] for i in range(n)]
        best = ranks[int(np.argmax(diff))]

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(ranks, errs, marker="o", label="Reconstruction error")
        ax.plot(ranks, interpol, alpha=0.7, marker="x", ls="--", label="Interpolated")
        ax.plot(ranks, diff, alpha=0.7, marker="^", ls="--", label="Difference")
        ax.axvline(best, color="red", alpha=0.25, ls="--", label=f"Elbow (rank {best})")

        ax.set_xlabel("Protein Rank")
        ax.set_ylabel("Reconstruction Error")
        ax.set_title("Reconstruction Error vs. Protein Rank")
        ax.set_xticks([k for k in ranks if k % 4 == 0 or k == 1])
        ax.legend()
        fig.tight_layout()
        return fig, ax

    def factor_correlation(self):
        self._ensure_global_phate()
        corr = np.corrcoef(self._global_scores.T)

        fig, ax = plt.subplots(figsize=(max(4, self.n_factors * 1.5), max(4, self.n_factors * 1.2)))
        vmax = np.abs(corr).max() or 1.0
        im = ax.imshow(corr, cmap=CMAP, vmin=-vmax, vmax=vmax, aspect="auto")

        ax.set_xticks(np.arange(self.n_factors))
        ax.set_xticklabels(self.factor_names, fontsize=9)
        ax.set_yticks(np.arange(self.n_factors))
        ax.set_yticklabels(self.factor_names, fontsize=9)

        self._annotate_heatmap(ax, corr, vmax)

        fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.5, pad=0.08)
        ax.set_title("Factor correlation (cell-level Pearson r)", fontsize=13)
        fig.tight_layout()
        return fig, ax

    # ------------------------------------------------------------------ #
    # PHATE-based visualizations
    # ------------------------------------------------------------------ #

    def phate_categorical(self):
        self._ensure_global_phate()

        phate_df = pd.DataFrame(
            self._global_phate, columns=["PHATE 1", "PHATE 2"],
            index=self._adata_sub.obs_names,
        )
        phate_df["treatment"] = self._adata_sub.obs[self.cct.treatment_col].astype(str).values
        phate_df["phase"] = self._adata_sub.obs[self.cct.cc_phase_col].astype(str).values

        dominant = self._adata_sub.to_df().idxmax(axis=1)
        phate_df["dominant_protein"] = dominant.str.removeprefix("nuclear_mean_").values

        configs = [
            ("treatment", "Treatment", {}),
            ("phase", "Cell-cycle phase", {}),
            ("dominant_protein", "Dominant protein", {"palette": "tab20"}),
        ]

        fig, axes = plt.subplots(1, 3, figsize=(16, 4))
        sns.despine(fig=fig, top=True, right=True, left=True, bottom=True)

        for ax, (col, title, kw) in zip(axes, configs):
            sns.scatterplot(
                data=phate_df, x="PHATE 1", y="PHATE 2", hue=col,
                s=self.dot_size, linewidth=0, alpha=0.8, ax=ax, **kw,
            )
            ax.set_title(title)
            ax.legend(title=title, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
            ax.set_xticks([])
            ax.set_yticks([])

        axes[0].set_xlabel("PHATE 1")
        axes[0].set_ylabel("PHATE 2")
        for ax in axes[1:]:
            ax.set_xlabel("")
            ax.set_ylabel("")

        fig.tight_layout()
        return fig, axes

    def phate_factor_scores(self):
        self._ensure_global_phate()

        rows, cols = self._grid_shape(self.n_factors, max_cols=2)
        fig, axes = plt.subplots(rows, cols, figsize=(8 * cols, 6 * rows))
        axes = np.atleast_1d(axes).ravel()
        sns.despine(fig=fig, top=True, right=True, left=True, bottom=True)

        for k, ax in enumerate(axes):
            if k >= self.n_factors:
                ax.axis("off")
                continue

            c = self._global_scores[:, k]
            vmax = np.nanpercentile(np.abs(c), 99) or 1.0

            pts = ax.scatter(
                self._global_phate[:, 0], self._global_phate[:, 1],
                c=c, cmap=CMAP, vmin=-vmax, vmax=vmax,
                s=self.dot_size, linewidth=0, alpha=0.8,
            )
            ax.set_title(self.factor_names[k])
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.set_xticks([])
            ax.set_yticks([])
            fig.colorbar(pts, ax=ax, shrink=0.7, pad=0.02)

        fig.tight_layout()
        return fig, axes

    def factor_score_violins(self):
        self._ensure_per_treatment()

        rows = []
        for t in self.treatments:
            pt = self._per_treatment[t]
            for f, fname in enumerate(self.factor_names):
                rows.append(pd.DataFrame({
                    "treatment": t,
                    "phase": pt["phase"],
                    "factor": fname,
                    "score": pt["embed"][:, f],
                }))
        df = pd.concat(rows, ignore_index=True)

        fig, axes = plt.subplots(
            self.n_factors, 1, figsize=(12, 2.4 * self.n_factors),
            sharex=True, squeeze=False,
        )
        axes = axes.ravel()

        for ax, fname in zip(axes, self.factor_names):
            sns.violinplot(
                data=df[df["factor"] == fname],
                x="treatment", y="score", hue="phase",
                order=list(self.treatments), hue_order=list(self.cc_phases),
                cut=0, inner="quartile", linewidth=0.8, ax=ax,
            )
            ax.set_title(fname, loc="left", fontsize=10)
            ax.axhline(0, color="0.6", lw=0.8, ls="--", zorder=0)
            ax.set_xlabel("")
            ax.set_ylabel("Factor score")
            if ax.get_legend() is not None:
                ax.get_legend().remove()
            sns.despine(ax=ax)

        axes[-1].tick_params(axis="x", rotation=45)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, title="Phase", bbox_to_anchor=(0.95, 0.5), loc="center left")
        fig.suptitle("Factor scores by treatment", y=0.96)
        fig.tight_layout(rect=[0, 0, 0.93, 0.99])
        return fig, axes

    def cosine_distance_violins(self):
        self._ensure_per_treatment()

        rows = []
        for t in self.treatments:
            pt = self._per_treatment[t]
            for f, fname in enumerate(self.factor_names):
                rows.append(pd.DataFrame({
                    "treatment": t,
                    "phase": pt["phase"],
                    "factor": fname,
                    "cos_dist": pt["cos_dist"][:, f],
                }))
        df = pd.concat(rows, ignore_index=True)

        fig, axes = plt.subplots(
            self.n_factors, 1, figsize=(12, 2.4 * self.n_factors),
            sharex=True, sharey=True, squeeze=False,
        )
        axes = axes.ravel()

        for ax, fname in zip(axes, self.factor_names):
            sns.violinplot(
                data=df[df["factor"] == fname],
                x="treatment", y="cos_dist", hue="phase",
                order=list(self.treatments), hue_order=list(self.cc_phases),
                cut=0, inner="quartile", linewidth=0.8, ax=ax,
            )
            ax.set_title(fname, loc="left", fontsize=10)
            ax.axhline(1.0, color="0.6", lw=0.8, ls="--", zorder=0)
            ax.set_xlabel("")
            ax.set_ylabel("Cosine distance")
            if ax.get_legend() is not None:
                ax.get_legend().remove()
            sns.despine(ax=ax)

        axes[-1].tick_params(axis="x", rotation=45)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, title="Phase", bbox_to_anchor=(0.95, 0.5), loc="center left")
        fig.suptitle("Cosine distance to factor archetype", y=0.96)
        fig.tight_layout(rect=[0, 0, 0.93, 0.99])
        return fig, axes

    def phate_local_by_factor(self, show_phase_contours=False):
        self._ensure_per_treatment()

        all_embeds = np.concatenate([self._per_treatment[t]["embed"] for t in self.treatments])
        vmax = np.nanpercentile(np.abs(all_embeds), 99) or 1.0
        nt = len(self.treatments)

        fig, axes = plt.subplots(
            self.n_factors, nt,
            figsize=(2.5 * nt, 2.2 * self.n_factors),
            squeeze=False,
        )

        for f in range(self.n_factors):
            for a, t in enumerate(self.treatments):
                ax = axes[f][a]
                pt = self._per_treatment[t]
                coords = pt["phate_coords"]
                pts = ax.scatter(
                    coords[:, 0], coords[:, 1],
                    c=pt["embed"][:, f], cmap=CMAP,
                    vmin=-vmax, vmax=vmax, s=self.dot_size, linewidth=0, alpha=0.8,
                )
                if show_phase_contours:
                    self._add_phase_contours(ax, coords, pt["phase"])
                ax.set_xticks([])
                ax.set_yticks([])
                if f == 0:
                    ax.set_title(t, fontsize=10)
                if a == 0:
                    ax.set_ylabel(self.factor_names[f], fontsize=10)

        title = "Per-treatment PHATE by latent protein-factor"
        if show_phase_contours:
            handles = [Line2D([], [], color=c, lw=1.5, label=p) for p, c in PHASE_COLORS.items()]
            fig.suptitle(title, y=0.98, x=0.38)
            fig.tight_layout(rect=[0, 0, 0.85, 0.99])
        else:
            fig.suptitle(title, y=0.98, x=0.38)


        cbar = fig.colorbar(pts, ax=axes, shrink=0.4, aspect=40, pad=0.02, label="Factor Value")
        if show_phase_contours:
            fig.legend(
                handles=handles, title="Phase", bbox_to_anchor=(cbar.ax.get_position().x0, 0.15),
                loc="upper left", fontsize=8, framealpha=0.7
            )

        return fig, axes

    def phate_global_by_factor(self, show_phase_contours=False):
        self._ensure_global_phate()

        _, embed_global = self._embed_factors(self._adata_sub.X)
        treatment_labels = self._adata_sub.obs[self.cct.treatment_col].astype(str).values
        phase_labels = self._adata_sub.obs[self.cct.cc_phase_col].astype(str).values
        vmax = np.nanpercentile(np.abs(embed_global), 99) or 1.0
        nt = len(self.treatments)

        fig, axes = plt.subplots(
            self.n_factors, nt,
            figsize=(2.5 * nt, 2.2 * self.n_factors),
            squeeze=False,
        )

        for f in range(self.n_factors):
            for a, t in enumerate(self.treatments):
                ax = axes[f][a]
                mask = treatment_labels == t
                coords = self._global_phate[mask]
                pts = ax.scatter(
                    coords[:, 0], coords[:, 1],
                    c=embed_global[mask, f], cmap=CMAP,
                    vmin=-vmax, vmax=vmax, s=self.dot_size, linewidth=0, alpha=0.9,
                )
                if show_phase_contours:
                    self._add_phase_contours(ax, coords, phase_labels[mask])
                ax.set_xticks([])
                ax.set_yticks([])
                if f == 0:
                    ax.set_title(t, fontsize=10)
                if a == 0:
                    ax.set_ylabel(self.factor_names[f], fontsize=10)

        title = "Per-treatment global PHATE embedding by latent protein-factor"
        if show_phase_contours:
            handles = [Line2D([], [], color=c, lw=1.5, label=p) for p, c in PHASE_COLORS.items()]
            fig.suptitle(title, y=0.98, x=0.38)
            fig.tight_layout(rect=[0, 0, 0.85, 0.99])
        else:
            fig.suptitle(title, y=0.98, x=0.38)

        cbar = fig.colorbar(pts, ax=axes, shrink=0.4, aspect=40, pad=0.02, label="Factor Value")
        if show_phase_contours:
            fig.legend(
                handles=handles, title="Phase", bbox_to_anchor=(cbar.ax.get_position().x0, 0.15),
                loc="upper left", fontsize=8, framealpha=0.7
            )

        return fig, axes

    # ------------------------------------------------------------------ #
    # Pseudotime visualizations
    # ------------------------------------------------------------------ #

    def pseudotime_factor_curves(self, factor=0, frac=0.15):
        self._ensure_per_treatment_dpt()

        nt = len(self.treatments)
        fname = self.factor_names[factor]

        row_h = 2.5
        fig, axes = plt.subplots(
            nt, 2, figsize=(10, row_h * nt), squeeze=False,
            gridspec_kw={"width_ratios": [1.4, 1], "hspace": 0.15, "wspace": 0.05},
        )

        ptime_max = max(
            np.nanmax(self._per_treatment[t]["dpt_pseudotime"])
            for t in self.treatments
            if np.any(np.isfinite(self._per_treatment[t]["dpt_pseudotime"]))
        )

        for row, t in enumerate(self.treatments):
            pt = self._per_treatment[t]
            ptime = pt["dpt_pseudotime"]
            scores = pt["embed"][:, factor]
            phases = pt["phase"]
            valid = np.isfinite(ptime)

            ax_curve = axes[row][0]
            ax_phate = axes[row][1]

            ptime_v = ptime[valid]
            scores_v = scores[valid]
            phases_v = phases[valid]

            for phase in ["G1", "S", "G2"]:
                m = phases_v == phase
                if m.any():
                    ax_curve.scatter(
                        ptime_v[m], scores_v[m],
                        color=PHASE_COLORS.get(phase, "0.5"),
                        s=3, alpha=0.15, rasterized=True,
                    )

            smooth = lowess(scores_v, ptime_v, frac=frac, return_sorted=True)
            ax_curve.plot(smooth[:, 0], smooth[:, 1], color="k", lw=2)
            ax_curve.axhline(0, color="0.6", lw=0.8, ls="--", zorder=0)

            ax_curve.set_ylabel(t, fontsize=10)
            ax_curve.set_xlim(0, ptime_max * 1.02)
            sns.despine(ax=ax_curve)
            if row < nt - 1:
                ax_curve.set_xticklabels([])

            coords = pt["phate_coords"]
            pts = ax_phate.scatter(
                coords[valid, 0], coords[valid, 1],
                c=ptime_v, cmap="viridis", vmin=0, vmax=ptime_max,
                s=self.dot_size, linewidth=0, alpha=0.8,
            )
            ax_phate.set_aspect("equal", adjustable="datalim")
            ax_phate.set_xticks([])
            ax_phate.set_yticks([])
            sns.despine(ax=ax_phate, left=True, bottom=True)
            if row == 0:
                ax_phate.set_title("PHATE (colored by pseudotime)", fontsize=10)

        axes[-1][0].set_xlabel("Pseudotime (DPT)")
        axes[0][0].set_title(f"{fname} score vs. pseudotime", fontsize=10)

        handles = [Line2D([], [], marker="o", ls="", color=c, ms=5) for c in PHASE_COLORS.values()]
        axes[0][0].legend(handles, PHASE_COLORS.keys(), loc="upper right", fontsize=8, framealpha=0.7)

        fig.suptitle(f"{fname} along diffusion pseudotime (per treatment)", y=0.91, fontsize=13)
        fig.tight_layout(rect=[0, 0, 0.90, 0.99])
        fig.colorbar(
            pts, ax=axes[:, 1].tolist(),
            shrink=0.25, aspect=25, pad=0.08, label="Pseudotime",
        )
        return fig, axes
