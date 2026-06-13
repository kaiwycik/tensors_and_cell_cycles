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

Let $\mathcal{X} \in \mathbb{R}^{I \times J \times K}$,  then the Tucker decomposition is:
$$
\begin{split}
\mathcal{X} \approx \mathcal{G} \times_{1} A \times_{2} B \times_{3} C & = \sum\limits_{p=1}^{P}\sum\limits_{q=1}^{Q}\sum\limits_{r=1}^{R} g_{pqr} \cdot \mathbf{a}_{p} \circ \mathbf{b}_{q} \circ \mathbf{c}_{r} \\
x_{ijk} & \approx \sum\limits_{p=1}^{P}\sum\limits_{q=1}^{Q}\sum\limits_{r=1}^{R} g_{pqr} \cdot a_{ip} \cdot b_{jq} \cdot c_{kr}
\end{split}
$$
where $\mathbf{A} \in \mathbb{R}^{I \times P}$, $\mathbf{B} \in \mathbb{R}^{J \times Q}$, and $\mathbf{C} \in \mathbb{R}^{K \times R}$ are factor matrices, $\mathbf{a}_p$, $\mathbf{b}_q$, and $\mathbf{c}_r$ are the respective column-vectors, and $\mathcal{G} \in \mathbb{R}^{P \times Q \times R}$. What this means is that each entry of $\mathcal{X}$ is a "column-wise sum" of the columns making up $\mathbf{A}$, $\mathbf{B}$, and $\mathbf{C}$, scaled by $\mathcal{G}$.

The matricized form of this is:
$$
\begin{split}
\mathbf{X}_{(1)} & \approx \mathbf{A}\mathbf{G}_{(1)}(\mathbf{C} \otimes \mathbf{B})^{T} \\
\mathbf{X}_{(2)} & \approx \mathbf{B}\mathbf{G}_{(2)}(\mathbf{C} \otimes \mathbf{A})^{T} \\
\mathbf{X}_{(3)} & \approx \mathbf{C}\mathbf{G}_{(3)}(\mathbf{B} \otimes \mathbf{A})^{T}, \\
\end{split}
$$
where $\mathbf{X}_{(n)}$ denotes the $n$-mode matricization of $\mathcal{X}$. Thinking through the first version of these, we matricize $\mathcal{G} \in \mathbb{R}^{P \times Q \times R}$ into $\mathbf{G}_{(1)} \in \mathbb{R}^{P \times (Q \cdot R)}$, to which we apply the first transformation $\mathbf{A} \mathbf{G}_{(1)} \in \mathbb{R}^{I \times (Q \cdot R)}$. Considering $(\mathbf{C} \otimes \mathbf{B}) \in \mathbb{R}^{(J \cdot K) \times (Q \cdot R)}$, we then obtain $\mathbf{A}\mathbf{G}_{(1)}(\mathbf{C} \otimes \mathbf{B})^{T} \in \mathbb{R}^{I \times (J \cdot K)}$.

### The $n$-Rank:
We define the $n$-rank ($R_{n} \stackrel{\Delta}{=} \mathit{rank}_n(\mathcal{X})$) as the dimension of the vector space spanned by the mode-$n$ fibers of $\mathcal{X}$ – we say that $\mathcal{X}$ is a rank-$(R_{1}, R_{2}, R_{3})$ tensor (in the case of $3^{\mathit{rd}}$ order tensors). Trivially, $R_{1} \leq I, R_{2} \leq J, R_{3} \leq K$.

Indeed, if we construct a Tucker decomposition with $R_{n} < \mathit{rank}_{n}(\mathcal{X})$, then our reconstruction will necessarily be inexact. We may exactly reconstruct $\mathcal{X}$ if we compose a Tucker decomposition with ranks equal to those of $\mathcal{X}$.

### Computing the Tucker Decomposition:
There are multiple ways of computing the Tucker decomposition. The idea of this first variant is to find the components that best capture the variation in mode $n$, independent of other modes.

![[HOSVD_algorithm.png]]

