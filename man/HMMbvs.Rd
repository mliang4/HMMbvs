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
HMMbvs(datmat=NULL, beta=NULL, gamma=NULL, nobs=NULL, sp=NULL, 
       hmm=TRUE, iterations=50000, pi=1, v1=1, v2=1, pt=31,
       pe=21, a=1, b=9, hidden_int=25, labelswitching=0)
}
\arguments{
  \item{datmat}{N x (1+pt+pe+2) matrix of observed data, the columns are: id, transition covariates, emission covariates, observed states and time difference delta. N = sum(nobs).}
  \item{beta}{(2pt+2pe) x 1 vector of initial regression coefficients}
  \item{gamma}{(2pt+2pe) x 1 vector of initial selection parameters}
  \item{nobs}{m x 1 vector of number of observations for each individual}
  \item{sp}{m x 1 vector of first observation positions}
  \item{hmm}{Indicator of which model to run, TRUE=run HMM, FALSE=run MSM}
  \item{iterations}{number of MCMC iterations}
  \item{pi}{Initial probability}
  \item{v1}{Prior variance}
  \item{v2}{Proposal variance}
  \item{pt}{Number of transition covariates (including intercept)}
  \item{pe}{Number of emission covariates (including intercept)}
  \item{a}{parameter for beta prior for inclusion probability}
  \item{b}{parameter for beta prior for inclusion probability}
  \item{hidden_int}{Interval of saving hidden states}
  \item{labelswitching}{Methods to deal with labelswitching. 0=doing nothing, 1=reject within step proposal if int0 < int1, 
  2=Swap intercepts of beta0 and beta1 if int0 < int1, 3=Swap beta0 and beta1 if int0 < int1}
}
\value{
  \code{\link{HMMbvs}} returns a list of the MCMC sampling results including the posterior samples of the coefficients and
  inclusion parameters, as well as the samples of the hidden chains and the log likelihood for each iteration.
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
# simulation study in paper

data = sim_data()
hmmDatmat = data$DATA # dataset with biased observation

InitBeta = c(rep(c(1,rep(0,pt-1)), 2), rep(c(1,rep(0,pe-1)), 2))
InitGamma = as.numeric(InitBeta!=0) 
InitGamma[c(1,pt+1,2*pt+1,2*pt+pe+1)]=1
  
# get obs length for each individual
len = NULL
for (i in 1:(nrow(hmmDatmat)-1)){
  if (hmmDatmat[i,1]!=hmmDatmat[i+1,1]){
    len = c(len,i)
  }
}
  
startIndex = c(0,len) # index for firstobs for each subject
endIndex = c(len,dim(hmmDatmat)[1]) # index for last obs for each subject
subjLen = endIndex-startIndex # number of observation for each subject
  
output = HMMbvs(datmat=hmmDatmat, beta=InitBeta, gamma=InitGamma, 
nobs=subjLen, sp=startIndex)
}
