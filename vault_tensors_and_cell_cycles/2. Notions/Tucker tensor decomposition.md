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

![[tucker_decomposition_example.png]]

Let $\mathcal{x} \in \mathbb{R}^{I \times J \times K}$,  then the Tucker decomposition is:
$$
\begin{split}
\mathcal{X} \approx \mathcal{G} \times_{1} A \times_{2} B \times_{3} C & = \sum\limits_{p=1}^{P}\sum\limits_{q=1}^{Q}\sum\limits_{r=1}^{R} g_{pqr} \cdot \mathbf{a}_{p} \circ \mathbf{b}_{q} \circ \mathbf{c}_{r} \\
x_{ijk} & \approx \sum\limits_{p=1}^{P}\sum\limits_{q=1}^{Q}\sum\limits_{r=1}^{R} g_{pqr} \cdot a_{ip} \cdot b_{jq} \cdot c_{kr}
\end{split}
$$
where $\mathbf{A} \in \mathbb{R}^{I \times P}$, $\mathbf{B} \in \mathbb{R}^{J \times Q}$, and $\mathbf{C} \in \mathbb{R}^{K \times R}$ are factor matrices and $\mathbf{a}_p$, $\mathbf{b}_q$, and $\mathbf{c}_r$ are the respective column-vectors. What this means is that each entry of $\mathcal{X}$ is a "column-wise sum" of the columns making up $\mathbf{A}$, $\mathbf{B}$, and $\mathbf{C}$, scaled by $\mathcal{G}$.



