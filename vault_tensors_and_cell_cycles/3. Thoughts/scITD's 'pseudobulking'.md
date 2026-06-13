---
category: thinking-note
---
# scITD's 'pseudobulking'
[[@Mitchel.etal2025]] perform "pseudobulking" before tensor decomposition.
Pseudobulking appears to just be grouping the cells and then normalizing & scaling by grouping (e.g., cell-types). 

>[!warning] !
> - What's unclear is why pseudobulk at all rather than leaving the data in single-cell, unlabelled form.
> - If one does have labels, aggregating by said labels is a reasonable way to group cells. But we would like to generalize to a case where we do not necessarily have labels.
