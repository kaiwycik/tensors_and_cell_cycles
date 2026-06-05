---
category: thinking-note
---
# On the type of tensor decomposition used in scITD

[[@Mitchel.etal2025]] write: "While other tensor decompositions can provide a similarly structured output, we chose the Tucker method for two main reasons. First, Tucker allows one to rotate its factors using approaches such as independent component analysis to improve interpretability of the results (Methods). Second, the raw output structure of Tucker gives the flexibility for the same sets of cell types to take part in separate multicellular patterns if there exist independent sets of covarying genes in those cell types. This would not be possible with other decompositions such as rank decomposition (that is, canonical polyadic decomposition), which would only produce one set of genes per set of cell types."

> [!hint] ?
> Are there other types of tensor decompositions which may provide different meaningful insight into interindividual patterns?
