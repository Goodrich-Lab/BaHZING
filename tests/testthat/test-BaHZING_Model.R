test_that("test BaHZING_Model", {
  ### Load example data
  data("iHMP_Reduced")

  ### Format microbiome data
  formatted_data <- Format_BaHZING(iHMP_Reduced)

  ### Perform Bayesian hierarchical zero-inflated negative binomial regression with g-computation
  ### Specify a mixture of exposures
  x <- c("soft_drinks_dietnum","diet_soft_drinks_dietnum"
         # "fruit_juice_dietnum","water_dietnum",
         # "alcohol_dietnum","yogurt_dietnum","dairy_dietnum","probiotic_dietnum","fruits_no_juice_dietnum",
         # "vegetables_dietnum","beans_soy_dietnum","whole_grains_dietnum","starch_dietnum"
         # ,"eggs_dietnum",
         # "processed_meat_dietnum","red_meat_dietnum","white_meat_dietnum","shellfish_dietnum","fish_dietnum",
         # "sweets_dietnum"
         )

  ### Specify a set of covariates
  # covar <- c("consent_age","sex",paste0("race",0:3),paste0("educ",0:7))
  covar <- c("consent_age")

  ###Set exposure standaridization to standard_normal or quantiles
  exposure_standardization = "standard_normal"

  ### Perform BaH-ZING with covariates
  results <- BaHZING_Model(formatted_data,
                           covar,
                           x,
                           exposure_standardization,
                           n.chains = 1,
                           n.adapt = 60,
                           n.iter.burnin = 2,
                           n.iter.sample= 2)
  # Test format data
  testthat::expect_equal(object = ncol(results), expected = 7)


  ### Perform BaH-ZING without covariates
  results <- BaHZING_Model(formatted_data = formatted_data,
                           covar = NULL,
                           x = x,
                           exposure_standardization = exposure_standardization,
                           n.chains = 1,
                           n.adapt = 60,
                           n.iter.burnin = 2,
                           n.iter.sample= 2)
  # Test format data
  testthat::expect_equal(object = ncol(results), expected = 7)
})
