---
category: notion-note
---
# Tucker tensor decomposition
Referring to [[@Kolda.Bader2009]], Tucker decompositions are a higher order generalization of SVD (or PCA).

> [!tldr] Summary
> Tucker decompositions reveal structure of interest about tensors, similar to how SVD reveals structure of interest about matrices.

> [!quote] Quote
> The Tucker decomposition is a form of higher-order PCA. It decomposes a tensor into a core tensor multiplied (or transformed) by a matrix along each mode.

_(For the sake of simplicity, we will only consider $3^{\text{rd}}$ order tensors here.)_

If we consider a tensor $\mathcal{X} \in \mathbb{R}^{I_1 \times I_2 \times I_3}$, its $n$-mode product with a matrix $U \in \mathbb{R}^{I_n \times J}$ is defined as:
$$

$$
