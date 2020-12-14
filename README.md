# Bayesian Variable Selection on HMM/MSM

This is an Rcpp package that performs Bayesian variable selection for mobile health data using Multistate Markov Model (MSM) or Hidden Markov Model (HMM). The main goal of this software package is to identify risk factors associated with individual behavior changes. The main function is written in C++ via the R package Rcpp, which allows both accelerated performance and easy implementation in R environment. This package also comes with R wrapper functions that produce processed outcome, summarization and plots for posterior inference. For more details, please see our manuscript

> Mingrui Liang, Matthew D. Koslovsky, Emily T. Hebert, Darla E. Kendzor, Michael S. Businelle, Marina Vannucci - *Bayesian Variable Selection for Binary Longitudinal Data with Measurement Error: An Application to mHealth Data*

## Getting Started

To get started, first make sure the following libraries are installed in R:

* Rcpp

* RcppArmadillo

* msm

* MASS

* coda

* reshape

* abind

* devtools

* ggplot2

* bayesplot


Switch to the desired directory to put the package, then copy and paste the following commands in terminal, which will download the package and switch to the working directory for you:

```shell
git clone https://github.com/mliang4/HMMbvs.git
cd HMMbvs
```

If the https repository doesn't work, the first line can be changed to:

```shell
git clone git@github.com:mliang4/HMMbvs.git
```

Then, open an R or RStudio console and set the working directory to the path of the package. Run the following command, which will install the package to the R library:

```R
library(devtools)
devtools::build(vignettes = F)
devtools::install()
```

## Using the package

For example of how to use this software package, please see vignettes/HMMbvs-vignette.html or the tutorial section in our paper.
