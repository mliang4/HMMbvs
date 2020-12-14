\name{HMMbvs-package}
\alias{HMMbvs-package}
\alias{HMMbvs}
\docType{package}
\title{
\packageTitle{HMMbvs}
}
\description{
\packageDescription{HMMbvs}
}
\usage{
HMMbvs_R( data = NULL, tcova = NULL, tforce = NULL, ecova = NULL, 
eforce = NULL, standardize = NULL, model = "HMM", init = "baseline", 
initvalue = NULL, iter = 5000, v1 = 5, v2 = 1, a = 1, b = 9, 
thin = 10, thin_hidden = 10)
}
\arguments{
  \item{datmat}{A long format matrix of observed data. Must contain at least three columns representing the subject ID, observed state, and observation time for each observation.}
  \item{tcova}{A named list or a vector representing covariates that may affect one of both directions of the transition}
  \item{ecova}{A named list or a vector representing covariates that may affect one of both directions of the emission}
  \item{tforce}{A named list representing covariates that are forced into the model}
  \item{eforce}{A named list representing covariates that are forced into the model}
  \item{standarize}{A vector of variable names to standardize}
  \item{model}{Either "HMM" or "MSM" indicating which model to run (Default = "HMM")}
  \item{init}{"baseline", "warmstart" or "manual", indicating MCMC initialization method (Default = "baseline")}
  \item{initvalue}{A named list indicating the value for each regression term to initialize}
  \item{iter}{Number of MCMC iterations (Default = 5000)}
  \item{v1}{Prior variance of the regression coefficients (Default = 5)}
  \item{v2}{Proposal variance of the regression coefficients (Default = 1)}
  \item{a}{Hyperparameter in beta prior for inclusion probability (Default = 1)}
  \item{b}{Hyperparameter in beta prior for inclusion probability (Default = 9)}
  \item{thin}{The number of MCMC samples are thinned to (Default = 10)}
  \item{thin_hidden}{The number of further thinning the simulated hidden states. For example, thin = 10, thin_hidden = 10 indicates the simulated hidden states are save every 100 iterations. (Default = 10)}
}
\value{
  \code{\link{HMMbvs}} returns a list including the MCMC samples of the coefficients and inclusion parameters, samples of the hidden chains, and the log-likelihood for each iteration.
}
\author{
\packageAuthor{HMMbvs}

Maintainer: \packageMaintainer{HMMbvs}
}
\references{
Mingrui Liang, Matthew D. Koslovsky, Emily T. Hebert, Darla E. Kendzor, Michael S. Businelle, Marina Vannucci - Bayesian Variable Selection for Binary Longitudinal Data with Measurement Error: An Application to mHealth Data
}

\keyword{ package }
\examples{
# Generating simulation data

set.seed(125)

df <- sim_data( ind = 150, window = 30, lambda_0 = 0.5, mu_0 = 0.5, 
emisP0 = 0.95, emisP1 = 0.95, pt = 30, pe = 20, bval = c(0.7,1.0,1.5), 
biased = TRUE )

# Sampling for the HMM model

hmmout <- HMMbvs_R( data = df$DATA, tcova = c( "tVar1", "tVar2", "tVar3", 
"tVar4", "tVar5", "tVar6", "tVar7", "tVar8", "tVar9", "tVar10", "tVar11", 
"tVar12", "tVar13", "tVar14", "tVar15" ), ecova = c( "eVar1", "eVar2", 
"eVar3", "eVar4", "eVar5", "eVar6", "eVar7" ), tforce = NULL, eforce = NULL, 
standardize = NULL, model = "HMM", init = "baseline", initvalue = NULL, 
iter = 25000, v1 = 5, v2 = 1, a = 1, b = 9, thin = 10, thin_hidden = 10)
  
# Accessing convergence and goodness-of-fit

convergence( output = hmmout, file = "convergence.pdf")

ypost = ppc( output = hmmout, file = "ppc_plot.pdf", type = "sum_one", 
burnin = 1000, postsample = NULL)

# Posterior inference of selection and estimation

out <- selection( output = hmmout, burnin = 1000, cred = 0.95, 
threshold = 0.5, trueGamma = NULL, plotting = TRUE, file = "mppi.pdf")
}
