---
category: notion-note
---
# scITD 
scITD is introduced in [[@Mitchel.etal2025|Coordinated, multicellular patterns of transcriptional variation that stratify patient cohorts are revealed by tensor decomposition (Mitchel et al., 2025)]].

> [!tldr] Summary
> The goal of scITD is to infer ‘multicellular patterns’ of gene expression. These multicellular patterns are sets of genes which covary together across individuals or donors in multiple cell types (not just in a single cell type).

> [!quote] Quote
> "Here, we developed an unsupervised method for analysis of interindividual variation, called single-cell interpretable tensor decomposition (scITD), that infers ‘multicellular patterns’ of gene expression in scRNA-seq datasets containing many source donors [...] We define a multicellular pattern to be a collection of genes in various cell types that covary together across donors. As our approach is unsupervised, it extracts the multicellular patterns with highest variance across donors and they may or may not be associated with outward donor phenotypes." 

To do so, the data is represented as a $n \times G \times C$ tensor, with $n \in \mathbb{N}$ being the number of donors, $G \in \mathbb{N}$ being the number of genes, and $C \in \mathbb{N}$ being the number of cell-clusters from pseudobulking.

Then: 
1. Within each cell type, but across individuals, the values are scaled and centered
2. A Tucker tensor decomposition is applied
3. The top $K$ informative factors are extracted. 

Each of the factors represents a multicellular gene expression pattern.