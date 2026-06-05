---
category: thinking-note
---
# scITD's 'pseudobulking'
[[@Mitchel.etal2025]] perform "pseudobulking" before tensor decomposition.
Pseudobulking appears to just be clustering the cells and then normalizing & scaling by cluster. The clusters act as a proxy for cell-type labels — so this only makes sense when one has (or is willing to infer) such labels.

>[!warning] !
> - What's unclear is why pseudobulk at all rather than leaving the data in single-cell, unlabelled form.
> - If one does have labels, aggregating by said labels is a reasonable way to group cells. But we would like to generalize to a case where we do not necessarily have labels.
