---
category: notion-note
---
# Tucker tensor decomposition
Referring to [[@Kolda.Bader2009]].

> [!tldr] Summary
> Tensors can be multiplied together, though obviously the notation and symbols for this are much more complex than for matrices.

> [!quote] Quote
> Here we consider only the tensor n-mode product, i.e., multiplying a tensor by a matrix (or a vector) in mode n.

_(For the sake of simplicity, we will only consider $3^{\text{rd}}$ order tensors here.)_

If we consider a tensor $\mathcal{X} \in \mathbb{R}^{I \times J \times K}$, its $2$-mode product with a matrix  $\mathbf{U} \in \mathbb{R}^{L \times J}$ is denoted as $\mathcal{Y} = \mathcal{X} \times_{n} \mathbf{U}$ and is of size $I \times L \times K$ (the $1$-mode and $3$-mode products are handled analogously).

In this case, for a fixed $i \in \{1, \dots, I\}$ and $k \in \{1, \dots, K\}$, we would obtain a fiber which "runs" along the columns of $\mathcal{X}$ – for $i=1, k=1$, this would be the front-facing top "edge" of our "cube".  Denote these fibers as $\mathbf{x}_{ik} \in \mathbb{R}^J$. 
We can think of $\mathcal{X} \times_{2} \mathbf{U}$ as multiplying all these fibers, for all possible combinations of $i, k$, with $\mathbf{U}$, yielding new fibers $\mathbf{y}_{ik} = \mathbf{U} \mathbf{x}_{ik} \in \mathbb{R}^L$. 

Analogously, we can think of "rotating" the individual fibers s.t. they are column vectors, and stacking all $\mathbf{x}_{ik}$ into a long matrix: $\mathbf{X} \in \mathbb{R}^{J \times (I \cdot K)}$ (see "Matricization", section 2.4 in [[@Kolda.Bader2009]]). We then obtain $\mathbf{Y} = \mathbf{U} \mathbf{X} \in \mathbb{R}^{L \times (I \cdot K)}$, which we can then again unfold into a tensor. 
This also provides some insight into the actual underlying computation:
$$
\begin{split}
\mathbf{Y}_{lm} & = \sum\limits_{j=1}^{J} u_{lj} \cdot x_{jm} \\
				& = \langle \mathbf{u}_{l}, \mathbf{x}_{m} \rangle,
\end{split}
$$
with $\mathbf{u}_l$ denoting the $l^{\text{th}}$ row of $\mathbf{U}$, and $\mathbf{x}_m$ denoting the $m^{\text{th}}$ column of $\mathbf{X}$, which is also the $m^{\text{th}}$ (in this case, row-wise) fiber of $\mathcal{X}$, with $m \in \{1, \dots, I \cdot J\}$. Or, analogously:
$$
\begin{split}
y_{ilk}   & = \sum\limits_{j = 1}^{J} u_{lj} \cdot x_{ijk} \\
					& = \langle \mathbf{u}_{l}, \mathbf{x}_{ik} \rangle
\end{split}
$$
A third perspective is to think of this as performing $k \in \{1, \dots, K \}$ matrix multiplications of $\mathbf{X}^{J \times I}$ slices of $\mathcal{X}$ which are subsequently stacked. 
_(Note that for this $2$-mode product, this can be thought of as "turning" over our "cube" of number onto its side, leaving the former columns now as rows and then considering the frontal slices.)_
We may index these slices as $\mathbf{X}_k$ and then considering:
$$
\mathbf{Y}_{k} = \mathbf{U} \mathbf{X}_k
$$


### Multiple $n$-mode products
If we again consider a tensor $\mathcal{G} \in \mathbb{R}^{I \times J \times K}$ and two matrices $\mathbf{U} \in \mathbb{R}^{L \times I}$ and $\mathbf{V} \in \mathbb{R}^{M \times J}$, we want to compute $\mathcal{X} = \mathcal{G} \times_{1} \mathbf{U} \times_{2} \mathbf{V}$. We can calculate $\mathcal{X}$ in two steps, first computing $\mathcal{Y} = \mathcal{G} \times_{1} \mathbf{U}$, and then computing $\mathcal{X} = \mathcal{Y} \times_{2} \mathbf{V}$. 

The first operation can be thought of by considering the frontal slices of $\mathcal{G}$, which we will denote as $\mathbf{G}_{k} \in \mathbb{R}^{I \times J}$, for $k \in \{ 1, \dots, K \}$. We then get:
$$
\mathbf{Y}_{k} = \mathbf{U} \mathbf{G}_{k} \in \mathbb{R}^{L \times J},
$$
which are now stacked into $\mathcal{Y}$. "Visually", the first coordinate (or $n^\text{th}$ coordinate in general) takes on the shape of $\mathbf{U}$'s first coordinate (in this case $L$) – to wit, our "cube" either "shrinks" or "grows" along the first (or $n^\text{th}$) coordinate. 
We then again consider the frontal slices of $\mathcal{Y}$, namely $\mathbf{Y}_{k} \in \mathbb{R}^{L \times J}$ – we will "rotate" (i.e., translate) these to be $\mathbf{Y'}_{k} \in \mathbb{R}^{J \times L}$ and can then again perform:
$$
\mathbf{X'}_{k} = \mathbf{V} \mathbf{Y'}_{k} \in \mathbb{R}^{M \times L},
$$
which provide us then with our $k$ slices $\mathbf{X}_k \in \mathbb{R}^{L \times M}$, which we can then stack to obtain our final tensor $\mathcal{X} \in \mathbb{R}^{L \times M \times K}$.
"Unfolding" what we did, each entry of $\mathcal{X}$ is given as 
$$
\begin{split}
x_{lmk}   & = \langle \mathbf{v}_{m}, \mathbf{y}_{lk} \rangle \\
					& = \sum\limits_{j = 1}^{J} v_{mj} \cdot y_{ljk} \\ 
					& = \sum\limits_{j = 1}^{J} v_{mj} \cdot \sum\limits_{i=1}^{I} u_{li} \cdot g_{ijk}
\end{split}
$$

We have now also worked our way up to understanding a multiplication of $\mathcal{G} \in \mathbb{R}^{P \times Q \times R}$ and three matrices $\mathbf{A} \in \mathbb{R}^{I \times P}$, $\mathbf{B} \in \mathbb{R}^{J \times Q}$, and $\mathbf{C} \in \mathbb{R}^{K \times R}$, namely:
$$
\begin{split}
x_{ijk} & = \sum\limits_{r=1}^{R} c_{kr} \cdot \sum\limits_{q = 1}^{Q} b_{jq} \cdot \sum\limits_{p=1}^{P} a_{ip} \cdot g_{pqr} \\
		& = \sum\limits_{r=1}^{R} \sum\limits_{q = 1}^{Q} \sum\limits_{p=1}^{P} a_{ip} \cdot b_{jq} \cdot c_{kr} \cdot g_{pqr} \\
\implies \mathcal{X} & = \sum\limits_{r=1}^{R} \sum\limits_{q = 1}^{Q} \sum\limits_{p=1}^{P} \mathbf{a}_{p} \circ \mathbf{b}_{q} \circ \mathbf{c}_{r} \circ g_{pqr},
\end{split}
$$
where $\mathbf{a}_p$, $\mathbf{b}_q$, and $\mathbf{c}_r$ are the respective column-vectors of $\mathbf{A}$, $\mathbf{B}$, and $\mathbf{C}$.