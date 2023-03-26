# Bayesian Multi-Study NMF for Mutational Signatures Analysis

We present a Bayesian multi-study NMF framework for mutational signatures analysis. Here, signatures can be discovered from multiple datasets, which might represent multiple studies, cancer types, conditions, or any other grouping structure of interest. In addition to the signatures themselves, we also explicitly model the sharing pattern of signatures across datasets, which allows each signature to potentially belong to any possible subset of the datasets. 

We developed two models in this framework: the *discovery-only* model, in which signatures are found in a fully unsupervised manner, and the *recovery-discovery* model, in which informative priors are used to encourage the recovery of previously known signatures, in addition to estimating novel signatures. We also developed *covariate extensions* of both models, which allows sample-level covariate effects to be simultaneously estimated alongside the signatures. 

Below, we describe the usage of our code. We also include example data and detailed demonstrations. 

## Usage: Discovery-Only Model

The only required input to the discovery-only model is a list of mutational counts matrices, where each such matrix corresponds to one dataset/group and is in *mutational motifs x samples* format. Optionally, there are hyperparameters that can be adjusted, but for most purposes, the default hyperparameters should suffice. If we call this list `M`, then the sampler can be run as follows:

```
source('Discovery_Only.R')
results <- discovery_sampler(M)
```

The output of this run is a list of three components. The first is a list of the most frequently occurring `A` matrices; the second is a list of corresponding signature estimates; and the third is a list of corresponding exposure estimates. The `A` matrix is a binary matrix indicating which signature belongs to each group, such that entry $(i,j) = 1$ if signature $j$ belongs to group $i$ and $0$ otherwise. To identify the best solution among these, the second step is to re-run the sampler conditional on each of these `A`s and estimate the marginal likelihood in each case. This can be done as follows:

```
results2 <- lapply(1:5,function(x)
  discovery_sampler(M,first=F,fixed=T,A.fixed=results[[1]][[x]],
                    inits=list(results[[2]][[x]],results[[3]][[x]])))
```

Note the additional arguments required here. The output of each individual call of `discovery_sampler` here is a list of three components: the estimated marginal log-likelihood, the fixed `A` matrix, and the posterior median estimate of the signatures matrix `P`. The solution with largest marginal log-likelihood can be thought of as the best solution. 

## Usage: Recovery-Discovery Model

The only required input to the recovery-discovery model is a list of mutational counts matrices, where each such matrix corresponds to one dataset/group and is in *mutational motifs x samples* format. Optionally, there are hyperparameters that can be adjusted, but for most purposes, the default hyperparameters should suffice. If we call this list `M`, then the sampler can be run as follows:

```
source('Recovery_Discovery.R')
results <- recovery_discovery(M)
```

The output of this run is a list of six components: the first three correspond to the discovery component (a list of the most frequently occurring `A` matrices in the discovery component, a list of corresponding signature estimates, and a list of corresponding exposure estimates) and the remaining three correspond to the recovery component (a list of the most frequently occurring `A` matrices in the recovery component, a list of corresponding signature estimates, and a list of corresponding exposure estimates). Here, the discovery component refers to novel detected signatures, and the recovery component refers to signatures with informative priors based on known signatures. The `A` matrix is a binary matrix indicating which signature belongs to each group, such that entry $(i,j) = 1$ if signature $j$ belongs to group $i$ and $0$ otherwise. In particular, the informative priors in this current version are built to represent the 78 signatures from the COSMIC v3.2 release. To identify the best solution among these, the second step is to re-run the sampler condition on these `A` matrices and estimate the marginal likelihood in each case. This can be done as follows:

```
results2 <- lapply(1:5,function(x)
  recovery_discovery(M,first=F,fixed=T,A.fixed=list(results[[1]][[x]],results[[4]][[x]]),
                     inits=list(results[[2]][[x]],results[[3]][[x]],results[[5]][[x]],results[[6]][[x]])))
```

Note the additional arguments required here. The output of each individual call of `recovery_discovery` here is a list of three components: the estimated marginal log-likelihood, a list of the fixed `A` matrices for the discovery and recovery component respectively, and a list of the posterior median estimates of the signatures matrix for the discovery and recovery component respectively. The solution with largest marginal log-likelihood can be thought of as the best solution. 

## Usage: Discovery-Only Model with Covariates

The two required inputs for the discovery-only model with covariates are a list of mutational counts matrices, where each such matrix corresponds to one dataset/group and is in *mutational motifs x samples* format, and a list of covariate matrices, where each such matrix corresponds to a dataset/group in the same order and is in *covariates x samples* format. Importantly, the covariates must be numeric! If you wish to encode any categorical covariates, this should be done using dummy variables. Optionally, there are hyperparameters that can be adjusted, but for most purposes, the default hyperparameters should suffice. If we call the list of mutational counts matrices as `M` and the list of covariate matrices as `cov`, then the sampler can be run as follows:
  
```
source('Discovery_Only_Covariates.R')
results <- discovery_covariates(M,cov)
```

The output of this run is a list of three components. The first is a list of the most frequently occurring `A` matrices; the second is a list of corresponding signature estimates; and the third is a list of corresponding exposure estimates. The `A` matrix is a binary matrix indicating which signature belongs to each group, such that entry $(i,j) = 1$ if signature $j$ belongs to group $i$ and $0$ otherwise. To identify the best solution among these, the second step is to re-run the sampler conditional on each of these `A`s and estimate the marginal likelihood in each case. This can be done as follows:

```
results2 <- lapply(1:5,function(x)
  discovery_covariates(M,cov,first=F,fixed=T,A.fixed=results[[1]][[x]],
                    inits=list(results[[2]][[x]],results[[3]][[x]])))
```

Note the additional arguments required here. The output of each individual call of `discovery_covariates` here is a list of five components: the estimated marginal log-likelihood, the fixed `A` matrix, the posterior median estimate of the signatures matrix `P`, a list of the posterior inclusion probabilities (PIPs) of each covariate for each signature in each dataset, and a list of the posterior median coefficient estimates of each covariate for each signature in each dataset. The solution with largest marginal log-likelihood can be thought of as the best solution. Note that the PIP and coefficient estimate matrices are in the format of *covariates x signatures*. Values of `NA` are used when a given signature does not belong to that dataset and, in the case of the coefficient estimate matrices, when a given covariate has a PIP of 0 for that signature.

## Usage: Recovery-Discovery Model with Covariates

The two required inputs for the recovery-discovery model with covariates are a list of mutational counts matrices, where each such matrix corresponds to one dataset/group and is in *mutational motifs x samples* format, and a list of covariate matrices, where each such matrix corresponds to a dataset/group in the same order and is in *covariates x samples* format. Importantly, the covariates must be numeric! If you wish to encode any categorical covariates, this should be done using dummy variables. Optionally, there are hyperparameters that can be adjusted, but for most purposes, the default hyperparameters should suffice. If we call the list of mutational counts matrices as `M` and the list of covariate matrices as `cov`, then the sampler can be run as follows:

```
source('Recovery_Discovery_Covariates.R')
results <- recovery_discovery_covariates(M,cov)
```

The output of this run is a list of six components: the first three correspond to the discovery component (a list of the most frequently occurring `A` matrices in the discovery component, a list of corresponding signature estimates, and a list of corresponding exposure estimates) and the remaining three correspond to the recovery component (a list of the most frequently occurring `A` matrices in the recovery component, a list of corresponding signature estimates, and a list of corresponding exposure estimates). Here, the discovery component refers to novel detected signatures, and the recovery component refers to signatures with informative priors based on known signatures. The `A` matrix is a binary matrix indicating which signature belongs to each group, such that entry $(i,j) = 1$ if signature $j$ belongs to group $i$ and $0$ otherwise. In particular, the informative priors in this current version are built to represent the 78 signatures from the COSMIC v3.2 release. To identify the best solution among these, the second step is to re-run the sampler condition on these `A` matrices and estimate the marginal likelihood in each case. This can be done as follows:

```
results2 <- lapply(1:5,function(x)
  recovery_discovery_covariates(M,cov,first=F,fixed=T,A.fixed=list(results[[1]][[x]],results[[4]][[x]]),
                    inits=list(results[[2]][[x]],results[[3]][[x]],results[[5]][[x]],results[[6]][[x]])))
```

Note the additional arguments required here. The output of each individual call of `recovery_discovery_covariates` here is a list of five components: the estimated marginal log-likelihood, a list of the fixed `A` matrices for the discovery and recovery component respectively, a list of the posterior median estimates of the signatures matrix for the discovery and recovery component respectively, a list of the lists of posterior inclusion probabilities (PIPs) of each covariate for each signature in each dataset in the discovery and recovery component respectively, and a list of the lists of posterior median coefficient estimates of each covariate for each signature in each dataset in the discovery and recovery component respectively. The solution with largest marginal log-likelihood can be thought of as the best solution. Note that the PIP and coefficient estimate matrices are in the format of *covariates x signatures*. Values of `NA` are used when a given signature does not belong to that dataset and, in the case of the coefficient estimate matrices, when a given covariate has a PIP of 0 for that signature.
