---
category: source-note
DOI: 10.1038/s41587-024-02411-z
title: Coordinated, multicellular patterns of transcriptional variation that stratify patient cohorts are revealed by tensor decomposition
authors: Jonathan Mitchel, M. Grace Gordon, Richard K. Perez, Evan Biederstedt, Raymund Bueno, Chun Jimmie Ye, Peter V. Kharchenko
year: 2025
itemType: journalArticle
publisher: Nature Publishing Group
tags: Gene regulatory networks, Gene expression
citekey: Mitchel.etal2025
annotation-target: file:///Users/kaiwycik/Zotero/storage/89IHZ67M/Mitchel%20et%20al.%20-%202025%20-%20Coordinated,%20multicellular%20patterns%20of%20transcriptional%20variation%20that%20stratify%20patient%20cohorts%20are%20r.pdf
status: read
dateread:
---
# Coordinated, multicellular patterns of transcriptional variation that stratify patient cohorts are revealed by tensor decomposition

> [!quote] Bibliography
> Mitchel, Jonathan, M. Grace Gordon, Richard K. Perez, et al. “Coordinated, Multicellular Patterns of Transcriptional Variation That Stratify Patient Cohorts Are Revealed by Tensor Decomposition.” _Nature Biotechnology_ 43, no. 7 (2025): 1192–201. [https://doi.org/10.1038/s41587-024-02411-z](https://doi.org/10.1038/s41587-024-02411-z).

> [!info]
> **Related**: 

> [!tldr] Summary
>
> Tissue-level and organism-level biological processes often involve the coordinated action of multiple distinct cell types. The recent application of single-cell assays to many individuals should enable the study of how donor-level variation in one cell type is linked to that in other cell types. Here we introduce a computational approach called single-cell interpretable tensor decomposition (scITD) to identify common axes of interindividual variation by considering joint expression variation across multiple cell types. scITD combines expression matrices from each cell type into a higher-order matrix and factorizes the result using the Tucker tensor decomposition. Applying scITD to single-cell RNA-sequencing data on 115 persons with lupus and 83 persons with coronavirus disease 2019, we identify patterns of coordinated cellular activity linked to disease severity and specific phenotypes, such as lupus nephritis. scITD results also implicate specific signaling pathways likely mediating coordination between cell types. Overall, scITD offers a tool for understanding the covariation of cell states across individuals, which can yield insights into the complex processes that define and stratify disease.
>.

> [!quote] Quotable
> .

## Notes
%% begin notes %%
**Introduction:**
- What is meant by “[…] principal components (PCs) from different cell types often show moderate, nonunique association.” ([Mitchel et al., 2025, p. 1193](zotero://select/library/items/VQRWTBPM))? Is the idea that the variation between cells across multiple individuals is not captured effectively by PCA?
- The goal of scITD is to infer ‘multicellular patterns’ of gene expression, i.e., patterns that don’t vary across cell types, but that vary across donors, individuals, &c. (“We define a multicellular pattern to be a collection of genes in various cell types that covary together across donors.” ([Mitchel et al., 2025, p. 1193](zotero://select/library/items/VQRWTBPM)))

**Results:**
- The fundamental idea is that the data is expressed as a tensor with shape $n \times G \times C$, where $n \in \mathbb{N}$ is the number of donors, $G \in \mathbb{N}$ is the number of genes, and $C \in \mathbb{N}$ is the number of cell-clusters from pseudobulking.
- Then: scale & center (in each cell type, across individuals) $\rightarrow$ Tucker tensor decomposition $\rightarrow$ get top $K$ informative factors (each representing a multicellular gene expression pattern).
- The authors then present a proof of concept with simulated data, two case studies, and biological interpretation of their results.
%% end notes %%


## Annotations
%% begin annotations %%
> we introduce a computational approach called single-cell interpretable tensor decomposition (scITD) to identify common axes of interindividual variation by considering joint expression variation across multiple cell types.

> Gene expression profiles derived from the same cell type vary across individuals, driven by a combination of genetics and environment.

> Furthermore, it has recently become more  common to generate single-cell RNA-sequencing (scRNA-seq) datasets from many individuals, allowing us to infer how the expression states of multiple cell types covary together across individuals.

> it can often be informative to identify the principal axes of variation stratifying individuals in an unsupervised

> When C cell populations from n donors are collapsed into

> (label-free) manner.

> it can be challenging to use such decompositions to identify coordinated variation between cells, as principal components (PCs) from different cell types often show moderate, nonunique association

> it can be highly informative to consider interindividual variation across multiple cell populations jointly.

> Here, we developed an unsupervised method for analysis of interindividual variation, called single-cell interpretable tensor decomposition (scITD), that infers ‘multicellular patterns’ of gene expression in scRNA-seq datasets containing many source donors (Fig. 1a).

> We define a multicellular pattern to be a collection of genes in various cell types that covary together across donors. As our approach is unsupervised, it extracts the multicellular patterns with highest variance across donors and they may or may not be associated with outward donor phenotypes.

> pseudobulk profiles, the dataset can be represented as a threedimensional matrix: a tensor T with dimensions n × G × C, where G is the number of genes (Fig. 1a, left).

> scale and center the expression of each gene (in each cell type) across individuals.

> apply the Tucker tensor decomposition

> extract the K most informative factors

> Each factor consists of a gene-by-cell type matrix of loading values, representing a multicellular gene expression pattern and a vector of sample scores, indicating the relative degree to which the respective pattern is expressed in each sample

> the sample score vector for a given factor

> describes each sample’s weighted average expression of genes in the corresponding loading matrix.

> While other tensor decompositions can provide a similarly structured output, we chose the Tucker method for two main reasons.

> First, Tucker allows one to rotate its factors

> Second, the raw output structure of Tucker gives the flexibility for the same sets of cell types to take part in separate multicellular patterns if there exist independent sets of covarying genes in those cell types.

> We introduced a computational approach called scITD, designed to identify common axes of interindividual gene expression variation by jointly considering gene expression states of multiple cell types. Such multicellular patterns can lend insight into mechanisms of tissue-level processes that are often relevant for disease. We principally applied our tool to large, clinically relevant, scRNA-seq datasets of PBMCs from persons affected by SLE or COVID-19, identifying multicellular patterns that stratified them.

%% end annotations %%


%% Import Date: 2026-06-05T11:49:14.454-04:00 %%
