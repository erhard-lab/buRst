# buRst <a href="https://github.com/erhard-lab/buRst"><img src="man/figures/logo.png" align="right" height="138" alt="buRst logo" /></a>

**Inferring transcriptional bursting dynamics from single-cell metabolic RNA-labeling data.**

buRst combines per-cell labeling likelihoods with mechanistic RNA count
distributions to infer burst frequency, burst size and RNA degradation rates.
It accounts for uncertainty in newly synthesized RNA measurements and supports
steady-state and transient models for perturbation experiments.

The package includes tools for:

- Single-gene and batch inference from scSLAM-seq data in `grandR` objects.
- Uncertainty estimation and joint comparisons of burst-frequency and burst-size changes.
- Simulation of RNA counts and metabolic-labeling evidence.
- Evaluation of theoretical RNA count distributions.

## Installation

Install the development version from GitHub. First install `rhdf5` through
[Bioconductor](https://bioconductor.org/packages/rhdf5/), then install buRst:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("rhdf5", quietly = TRUE)) {
  BiocManager::install("rhdf5", ask = FALSE, update = FALSE)
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("erhard-lab/buRst")
library(buRst)
```

Installing from source requires a suitable C++ toolchain. Other required R
packages are installed as dependencies. The package workflows have been tested
on Linux; memory requirements and runtime depend on the dataset and number of
genes. Sequential fitting is available with `parallel = FALSE`.

## Getting started

Choose the vignette that fits your next step.

| Goal | Vignette |
|---|---|
| See how different bursting models can have the same mean RNA abundance | [Same mean, different bursts](vignettes/bursting_same_mean.Rmd) |
| Generate labeled RNA counts and discrete likelihoods | [Simulate labeled RNA counts](vignettes/simulation.Rmd) |
| Follow a complete simulation and recover a known mechanism | [Find the mechanism of an acute perturbation](vignettes/acute_pertubation.Rmd) |
| Prepare your own GRAND3 inputs and interpret condition comparisons | [Fit bursting kinetics to scSLAM-seq data](vignettes/inference.Rmd) |

After installing with built vignettes, use `browseVignettes("buRst")` to open
rendered tutorials. GitHub installations omit vignette builds by default;
use `build_vignettes = TRUE` after installing `knitr` and `rmarkdown` to include them.

Function documentation is available from R:

```r
help(package = "buRst")
?run_kinetics_batch
?compare_posteriors_gate_star_kmB
?simulateBursting
```

Real-data inference uses processed GRAND3 count data and the corresponding
full discrete labeling-likelihood HDF5 file. Data preparation, experimental
conditions and labeling intervals are described in the inference vignette.

## Manuscript analyses

Analysis workflows, processed inputs, saved simulations, fitted results and
figures have reserved Zenodo DOI `10.5281/zenodo.24362854`. 
## Questions and bug reports

Please use the [GitHub issue tracker](https://github.com/erhard-lab/buRst/issues).
For reproducible bug reports, include a minimal example, the error message and
`sessionInfo()`.
