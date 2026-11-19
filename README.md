# Do bird eggs revert in shape after the loss of flight?

## Background

Stoddard et al. (2017, *Science*) found that birds that fly well tend to lay
more elongated, asymmetric eggs, and suggested that adaptations for flight
shaped the egg. Flight has been lost many times in birds (kiwi, kakapo,
rails, penguins and others). If flight really drives egg shape, those
flightless lineages should have drifted back towards rounder eggs. This
project tests that idea.

Short answer: I found no sign of such a reversion, although with only a
handful of independent flight losses a small effect cannot be ruled out.
The write-up is in [`report/report.pdf`](report/report.pdf).

## Methods

- **Data:** egg shape from Montgomerie et al. (2021; 946 species) and
  Stoddard et al. (2017; 1,388 species), analysed separately; flight status
  from Sayol et al. (2020); locomotion type from Saitta et al. (2025); body
  mass and hand-wing index from AVONET; phylogeny from BirdTree (100
  posterior trees and a consensus tree).
- **Phylogenetic regression (PGLS):** egg shape ~ flight status, controlling
  for body mass and egg size.
- **Sister-pair comparisons:** each flightless lineage against its closest
  flying relatives, compared with flying sister pairs of similar age.
- **Evolutionary rates:** multi-rate Brownian motion and OU models, checked
  against simulated data with no real rate difference.
- **Robustness:** leaving out one lineage at a time, alternative ratite
  trees, different body-mass models, and equivalence tests for how large
  an effect the data can exclude.

## Running it

`Rscript run_all.R` runs the whole pipeline in `R/` and renders the report.
Three data files have to be downloaded by hand first; `R/01_download.R`
says which ones and where to put them.

Code: MIT (`LICENSE`). Text and figures: CC BY 4.0.
