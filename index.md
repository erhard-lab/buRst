# buRst

### Inferring transcriptional bursting dynamics from single-cell metabolic RNA-labeling data.

**buRst** is an R package for likelihood-based inference of transcriptional bursting dynamics from single-cell RNA sequencing data with metabolic RNA labeling.

It combines mechanistic RNA count models with per-cell labeling likelihoods, allowing uncertainty in the measurement of newly synthesized RNA to be incorporated into the inference of transcriptional kinetic parameters.

## What can buRst do?

**Transcriptional bursting inference**

Estimate transcriptional burst frequency, burst size, and RNA degradation rates using mechanistic models of RNA production and degradation.

**Uncertainty-aware inference**

Account for uncertainty in newly synthesized RNA measurements through likelihood-based modeling of metabolic RNA labeling data.

**Comparison of experimental conditions**

Investigate changes in transcriptional kinetics between experimental conditions using posterior-based comparisons.

**Simulation**

Simulate RNA count distributions and labeling evidence under steady-state and transient transcriptional bursting models.

## Installation

Once buRst is available on CRAN, install it using:

```r
install.packages("buRst")
```

Until then, the development version can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("erhard-lab/buRst")
```

## Getting started

New to buRst? Start with the inference tutorial, which introduces the analysis workflow.

**[Getting started with inference](articles/inference.html)**

## Tutorials and vignettes

Explore the different aspects of transcriptional bursting analysis:

### Inference

Learn how to infer transcriptional bursting parameters from single-cell metabolic RNA labeling data.

[Read the inference vignette →](articles/inference.html)

### Simulation

Explore the simulation of RNA counts and labeling evidence under transcriptional bursting models.

[Read the simulation vignette →](articles/simulation.html)

### Same Mean, Different Bursts

Discover how different combinations of burst frequency and burst size can generate the same average RNA abundance while producing different RNA count distributions.

[Explore transcriptional bursting →](articles/bursting_same_mean.html)

### Acute Perturbation

Explore transcriptional bursting dynamics following an acute perturbation using simulated data.

[Read the acute perturbation vignette →](articles/acute_pertubation.html)

## Function reference

Looking for documentation of a particular function?

The complete function reference contains the available distributions, inference methods, simulation utilities, and data-processing functions.

**[Browse the function reference →](reference/index.html)**

## Source code

buRst is developed by the Erhard Lab.

The source code, development history, and issue tracker are available on GitHub:

[github.com/erhard-lab/buRst](https://github.com/erhard-lab/buRst)