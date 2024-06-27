# Bayesian Multi-Study NMF for Mutational Signatures Analysis

We present a Bayesian multi-study NMF framework for mutational signatures analysis. Here, signatures can be discovered from multiple datasets, which might represent multiple studies, cancer types, conditions, or any other grouping structure of interest. In addition to the signatures themselves, we also explicitly model the sharing pattern of signatures across datasets, which allows each signature to potentially belong to any possible subset of the datasets.

We developed two models in this framework: the *discovery-only* model, in which signatures are found in a fully unsupervised manner, and the *recovery-discovery* model, in which informative priors are used to encourage the recovery of previously known signatures, in addition to estimating novel signatures. We also developed *covariate extensions* of both models, which allows sample-level covariate effects to be simultaneously estimated alongside the signatures.

## Installation

The package requires the dependency `armspp`, and both can be installed as follows:

```
install.packages('armspp')
devtools::install_github("igrabski/MultiStudyNMF", build_vignettes = FALSE)
```

## Demos and Vignettes

We provide two basic datasets and show examples of their analysis in the vignettes (located in the "vignettes" folder). The datasets are from PCAWG, and can be accessed as follows:

```
data(lung)
data(lymph)
```

The lung dataset consists of two studies (lung adenocarcinoma and lung squamous cell carcinoma). The lymph dataset consists of two studies (representing two different cohorts) that are each a combination of Lymph BNHL and Lymph CLL; this dataset also comes with a covariate denoting the cancer ID. 

## Usage: Discovery-Only Model

The only required input to the discovery-only model is a list of mutational counts matrices, where each such matrix corresponds to one dataset/group and is in *mutational motifs x samples* format. It is not required, but can be helpful if the ordering of the mutational motifs matches that of the COSMIC reference data (see `data(cosmic)`) for downstream interpretation. If we call this list `M`, then in the first step, the sampler can be run as follows:

```
data(lung)
M <- lung$M
results <- discovery_fit(M)
```

If the mutational motifs of the input data matched the order of the COSMIC reference data (see `data(cosmic)`), we can use the `label_A` function to annotate the signature sharing matrix by its closest match:

```
label_A(results$A,results$P)
```

Optionally, we can also view the proportion of samples in each study with high exposure to each signature:

```
proportion_A(results,thresh=0.05,label=T) # set label = F if you don't want to annotate signatures
```

## Usage: Recovery-Discovery Model

The only required input to the recovery-discovery model is a list of mutational counts matrices, where each such matrix corresponds to one dataset/group and is in *mutational motifs x samples* format. *Important*: the order of mutational motifs has to match that of the COSMIC v3.2 reference, which can be viewed as `data(cosmic)`. If we call this list `M`, then the sampler can be run as follows:

```         
data(lung)
M <- lung$M
results <- recovery_discovery_fit(M)
```

We can again use the `label_A` and `proportion_A` helper functions (see vignettes for more detailed examples).

## Usage: Discovery-Only Model with Covariates

The two required inputs for the discovery-only model with covariates are a list of mutational counts matrices, where each such matrix corresponds to one dataset/group and is in *mutational motifs x samples* format, and a list of covariate matrices, where each such matrix corresponds to a dataset/group in the same order and is in *covariates x samples* format. Importantly, the covariates must be numeric! If you wish to encode any categorical covariates, this should be done using dummy variables. The same covariates in the same order should be supplied for each dataset/group.

If we call the list of mutational counts matrices as `M` and the list of covariate matrices as `cov`, then the sampler can be run as follows:

```         
data(lymph)
M <- lymph$M
cov <- lymph$cov
results <- discovery_cov_fit(M,cov)
```

We can again use the `label_A` and `proportion_A` helper functions (see vignettes for more detailed examples).

## Usage: Recovery-Discovery Model with Covariates

The two required inputs for the recovery-discovery model with covariates are a list of mutational counts matrices, where each such matrix corresponds to one dataset/group and is in *mutational motifs x samples* format, and a list of covariate matrices, where each such matrix corresponds to a dataset/group in the same order and is in *covariates x samples* format. *Important*: the order of mutational motifs has to match that of the COSMIC v3.2 reference, which can be viewed as `data(cosmic)`.

Additionally, the covariates must be numeric! If you wish to encode any categorical covariates, this should be done using dummy variables. The same covariates in the same order should be supplied for each dataset/group.

If we call the list of mutational counts matrices as `M` and the list of covariate matrices as `cov`, then the sampler can be run as follows:

```      
data(lymph)
M <- lymph$M
cov <- lymph$cov
results <- recovery_discovery_cov_fit(M,cov)
```

We can again use the `label_A` and `proportion_A` helper functions (see vignettes for more detailed examples).
