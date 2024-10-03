
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BaHZING: Bayesian Hierarchical Zero-Inflated Negative Binomial Regression with G-Computation

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/Goodrich-Lab/BaHZING/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Goodrich-Lab/BaHZING?branch=main)
[![R-CMD-check](https://github.com/Goodrich-Lab/BaHZING/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Goodrich-Lab/BaHZING/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Contents

- [Overview](#overview)
- [Installation](#overview)
- [Example](#example)

## Overview

This platform is dedicated to the analysis of environmental mixture
exposures and other mixture exposures with hierarchical gut microbiome
outcome data. BaHZING, short for Bayesian Hierarchical Zero-Inflated
Negative Binomial Regression with G-Computation, is a powerful toolkit
designed to uncover intricate relationships between environmental
exposures and gut microbiome outcomes.

## Installation

You can install the development version of BaHZING from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("Goodrich-Lab/BaHZING")
```

## Input Requirements

BaHZING’s formatting function requires a phyloseq object as the input.
Ensure that your phyloseq object adheres to the required format. Refer
to the phyloseq documentation for details on creating and manipulating
phyloseq objects. For covariates in phyloseq object: All categorical
covariates must be formatted as binary indicator variables rather than
categorical or factor variables. For missing data: Ensure that there are
no NA values in your data, perform the necessary filtering or imputation
methods necessary to handle these values as BaHZING cannot function with
missing data.

## Output

The formatting function produces an object that contains a table with
all exposure variables, all covariates and the lowest level of
microbiome outcome data. Additionally, the function produces binary
matrices that represent the relationship between each level of
microbiome data and the level above it.

The BaHZING model function provides results from the Bayesian
hierarchical zero-inflated negative binomial regression model.

## Example

For a quick start, we utilize data from the iHMP publicly available data
set. This example will guide you through the process of using BaHZING
for analyzing environmental mixture exposures with hierarchical gut
microbiome outcome data.

``` r
library(BaHZING)
#> Loading required package: rjags
#> Loading required package: coda
#> Linked to JAGS 4.3.2
#> Loaded modules: basemod,bugs

### Load example data
data("iHMP_Reduced")

### Format microbiome data
formatted_data <- Format_BaHZING(iHMP_Reduced)

### Perform Bayesian hierarchical zero-inflated negative binomial regression with g-computation
### Specify a mixture of exposures
x <- c("soft_drinks_dietnum","diet_soft_drinks_dietnum")
# "fruit_juice_dietnum","water_dietnum",
# "alcohol_dietnum","yogurt_dietnum","dairy_dietnum","probiotic_dietnum",
# "fruits_no_juice_dietnum", "vegetables_dietnum","beans_soy_dietnum", 
# "whole_grains_dietnum","starch_dietnum", "eggs_dietnum", 
# "processed_meat_dietnum","red_meat_dietnum","white_meat_dietnum",
# "shellfish_dietnum","fish_dietnum", "sweets_dietnum"

### Specify a set of covariates
covar <- c("consent_age")

###Set exposure standaridization to standard_normal or quantiles
exposure_standardization = "standard_normal"

### Perform BaH-ZING
# NOTE: For this example, we are using a small number of iterations to reduce 
# the runtime. For a real analysis, we recommend using n.chain = 3, 
# n.adapt, n.burnin, and n.sample > 5000.
results <- BaHZING_Model(formatted_data,
                         covar,
                         x,
                         exposure_standardization,
                         n.chains = 1,
                         n.adapt = 60,
                         n.burnin = 2,
                         n.sample = 2)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 23310
#>    Unobserved stochastic nodes: 26312
#>    Total graph size: 337340
#> 
#> Initializing model
```
