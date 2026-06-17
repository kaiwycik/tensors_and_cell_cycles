---
category: thinking-note
---
# scITD's 'pseudobulking'
[[@Mitchel.etal2025]] perform "pseudobulking" before tensor decomposition.
Pseudobulking appears to just be grouping the cells and then normalizing & scaling by grouping (e.g., cell-types). 

>[!warning] !
> - What's unclear is why pseudobulk at all rather than leaving the data in single-cell, unlabelled form.
> - If one does have labels, aggregating by said labels is a reasonable way to group cells. But we would like to generalize to a case where we do not necessarily have labels.

## Why not keep the samples as their own mode?
In [[2026-06-15 - Meeting w. Tarek|the 2026-06-15 meeting]] we considered avoiding pseudobulking by promoting the individual cells (samples) to a fourth tensor mode: $\mathcal{X} \approx \mathcal{G} \times_1 T \times_2 C \times_3 P \times_4 S$, of shape $\textit{treatment} \times \textit{phase} \times \textit{protein} \times \textit{sample}$. Having tried this out, we can now articulate why this is ill-posed.

A cell is measured exactly once. Its treatment and phase are *labels attached to it*, not axes it ranges along — each cell sits in exactly one $(\textit{treatment}, \textit{phase})$ combination. 

A proper $t \times c \times p \times s$ tensor, by contrast, assumes every sample exists under *every* $(\textit{treatment}, \textit{phase})$ pair. With $t = 8$ treatments and $c = 3$ phases (after excluding G0/M), each cell occupies only $1$ of $t \cdot c = 24$ slices. So the dense tensor would be:
- only $\tfrac{1}{t \cdot c} = \tfrac{1}{24} \approx 4.2\%$ filled — i.e. **~95.8% `NaN`**, "every now and then" a real entry;
- not even rectangular: the $(\textit{treatment}, \textit{phase})$ groups are unbalanced (Dinaciclib ≈ 23.5% of cells, Talazoparib ≈ 5%, …), so the $231{,}934$ cells don't divide evenly into $24$ groups ($231{,}934 / 24 = 9663.9\ldots$). One would have to `NaN`-pad the smaller groups up to the largest.

>[!fail] ...
> This is why the code throws a `reshape` error. `adata_sans_G0_M.X` holds $231{,}934 \times 13 = 3{,}015{,}142$ values; the requested shape $(8, 3, 13, 231{,}934)$ wants $72{,}363{,}408 = 24 \times$ as many. `reshape` cannot fabricate the missing 23/24. 

>[!hint] ?
> Pseudobulking is not optional for now. The sample axis is precisely the thing aggregation dissolves, because individual cells are not a shared, alignable unit across conditions: "cell #5000 under Dinaciclib/G1" has nothing to do with "cell #5000 under Control/S".

This sharpens the original question above: we may still want to generalize pseudobulking to the label-free case, but we cannot sidestep it by making samples a mode. A separate thought — [[Recovering sample-level resolution from Tucker factor loadings]] — asks whether sample-level resolution can instead be *recovered* after the fact from the decomposition.
