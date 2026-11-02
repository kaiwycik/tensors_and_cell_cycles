# Tensor analyses with the 96-feature set

This folder reruns every tensor analysis in `code/` on 96 features instead of 22. The data, preprocessing, decomposition and plots are otherwise the same.

## Feature set

The 96 features are:
- 49 nuclear: the 47 nuclear marker means (`nuclear_mean_<marker>`), `nuclear_area` and `DNA_int`
- 47 cytoplasmic: the 47 cytoplasmic marker means (`cyto_mean_<marker>`)

Excluded are NaKATPase, pRB and RB in both compartments, and the derived `pRB/RB`. Together with these 7, the 96 make up all 103 columns of each `<line>_merged.h5ad`. In the G1/S/G2 cells of every line, none of the 96 has missing values or zero variance, including within the Control cells used for Control referencing.

The 47 markers are listed in the `MARKERS` cell of each notebook. Features are always in the same order: nuclear markers, then `nuclear_area` and `DNA_int`, then cytoplasmic markers, with markers in alphabetical order within each block. Plots label them `nuc_<marker>` and `cyto_<marker>`, so the two measurements of one marker can be told apart. `nuclear_area` and `DNA_int` keep their names.

## What is the same as in `code/`

The pipeline is unchanged:
- pseudobulk means per group, z-scored per feature
- Control referencing in the per-drug notebooks
- HOSVD with full treatment (or cell-line) and phase ranks
- 4 latent protein-factors, as in `code/`, so the two analyses can be compared; the elbow rank from `determine_rank()` is still printed
- the same figures, the cross-cell-line sign alignment, and the CARVE clustering

Plot differences:
- The protein Gram clustermap has no cell values, because a 96 x 96 grid is too dense to annotate.
- The cross-cell-line loading violins (section 2.4) are split into four figures, two for the nuclear and two for the cytoplasmic features.
- In the per-cell-line PHATE overview, the dominant-protein panel colors the nine most frequent dominant features and groups the rest as "other".

The protein mode now has more features (96) than the tensor has pseudobulk groups (8 x 3 = 24). Its unfolding is therefore 96 x 24 and, after z-scoring, has at most 23 nonzero singular values. In the Control per-drug tensor the maximum is 16, because each line is referenced to its own Control cells. With 22 features this limit only applied to the Control tensor. Beyond the maximum, the reconstruction error in `determine_rank()` is zero, so the reported elbow mostly falls at the maximum and should not be read as a rank estimate.

The cross-cell-line notebook keeps two choices made from the 22-feature results. Sections 3.3 and 3.4 leave out Factors 1 and 4, and sections 4.1 and 4.2 show Factors 1 and 2. With 96 features the factors are different, so these choices may need revisiting.

Pseudotime builds its neighbor graph on the z-scored features (`sc.pp.neighbors(..., use_rep="X")` in `cc_tensor_viz.py`). Without this, scanpy would first reduce more than 50 features to 50 principal components. With 22 features scanpy already used the features directly, so the method is the same as before.

## Notebooks and run order

The notebooks import `cc_tensor` and `cc_tensor_viz` from `code/`, read the data from `code/data/`, and write results to `results/` in this folder.

1. Per cell line: `analysis_<line>.ipynb`. Each writes `results/cct_<line>.pkl`.
2. Per drug: `analysis_drug_<drug>.ipynb`. Each writes `results/cct_drug_<drug>.pkl`. These are independent of step 1.
3. `analysis__cross_cellline.ipynb`, which needs the eight pickles from step 1.
