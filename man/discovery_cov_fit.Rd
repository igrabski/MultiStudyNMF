% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Discovery_Only_Covariates.R
\name{discovery_cov_fit}
\alias{discovery_cov_fit}
\title{Run the Discovery-Only Sampler with Covariates}
\usage{
discovery_cov_fit(
  M,
  cov,
  lambda_p = 0.5,
  a_p = 5,
  b_p = 0.05,
  lambda_e = 0.1,
  a_e = 5,
  b_e = 0.01,
  Ntot = 50
)
}
\arguments{
\item{M}{List of matrices corresponding to each study, in features-by-samples format.}

\item{cov}{List of covariates matrices corresponding to each study,
in covariates-by-samples format; the same covariates in the same order should be supplied for
each study. Covariates must be numeric.}

\item{lambda_p}{Lambda_p, set to 0.5 by default.}

\item{a_p}{a_p, set to 5 by default.}

\item{b_p}{b_p, set to 0.05 by default.}

\item{lambda_e}{lambda_e, set to 0.1 by default.}

\item{a_e}{a_e, set to 5 by default.}

\item{b_e}{b_e, set to 0.01 by default.}

\item{Ntot}{Total possible number of signatures to consider, set to 50 by default.}
}
\value{
A named list describing the best solution:
\itemize{
 \item{BIC: }{the BIC value of the solution.}
 \item{A: }{the point estimate of the A matrix, which indicates whether each signature is
 present in each study.}
 \item{P: }{the posterior median of the signatures matrix.}
 \item{E: }{a list of the posterior medians of the exposures matrices for each study s,
 subsetted to just the signatures indicated by A to belong to that study.}
 \item{indics: }{a list of the posterior inclusion probabilities for each study s,
 describing the probability that each covariate affects each signature indicated by A
 to belong to that study.}
 \item{coefs: }{a list of the coefficients, conditional on inclusion, for each study s,
 describing the contribution of each covariate to each signature indicated by A to belong
 to that study.}
 \item{W: }{a list of the normalizing constants for each study; note that the input counts
 matrix M_s for each study s can be modeled as Poisson(W_s * P_s E_s), where P_s and E_s
 contain just the signatures that belong to study s.}
 }
}
\description{
This runs the fully unsupervised discovery-only sampler for multi-study NMF, with covariates.
}
