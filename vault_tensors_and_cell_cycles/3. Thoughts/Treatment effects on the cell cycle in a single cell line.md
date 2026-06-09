---
category: thinking-note
---
# Treatment effects on the cell cycle in a single cell line

Idea from [[2026-06-09|the meeting on 9th of June 2026]] – rather than reproducing [[scITD]] across donors, we consider a single cell line (start with T47D, possibly MCF7) and ask how *treatments* reshape the cell cycle. We may pseudobulk and scale/centre similarly to [[scITD's 'pseudobulking']], but consider the tensor as
$$
\textit{treatment} \times \textit{cell-cycle phase} \times \textit{protein}.
$$
A Tucker decomposition then gives
$$
X \approx G \times_1 T \times_2 C \times_3 P,
$$
and the object of interest are the loadings — to wit, $T$ (treatment) and $C$ (cell-cycle phases).

We want to answer three questions:
1. Which cell-cycle structure remains **unchanged** under treatment?
2. How are the cell-cycle phases **distributed** across treatments?
3. How do treatments **affect** the cell cycle?

> [!hint] ?
> What mode of variation does each treatment component ($T$) capture? A component that loads roughly evenly across treatments may point at cell-cycle structure that is invariant to treatment (question 1); one that splits the treatments points at a treatment-driven effect (question 3).

> [!faq] ?
> How do we quantify the effect of treatment components $T$ over cell-cycle phases? The coupling lives in the core $G$: a slice $G[\,:, c, p\,]$ says how strongly each treatment component participates for a given phase/protein pairing, so we can track how a $T$-component's weight moves across the $C$ and $P$ modes.

> [!warning] !
> One cell line first — so we cannot yet see how the loadings vary *across* cell lines (may expand in the future to a 4th mode). 
> Cautios for imbalance: treatments differ in cell counts, which may bias the pseudobulk means and hence the within-phase scaling.

Anchored in [[@Mitchel.etal2025]] (the [[scITD]] pipeline) and [[@Kolda.Bader2009]] for Tucker; data from [[@Zikry.etal2024]]. Relates to [[On the type of tensor decomposition used in scITD]].

