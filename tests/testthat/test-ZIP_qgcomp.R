test_that("test ZIP_qgcomp", {
  ### Load example data
  data("iHMP_Reduced")

  ### Format microbiome data
  formatted_data <- Format_BaHZING(iHMP_Reduced)

  ### Perform qgcomp with zero-inflated negative binomial regression/poisson regression
  ### Specify a mixture of exposures
  x <- c("soft_drinks_dietnum",
         "diet_soft_drinks_dietnum"
         # "fruit_juice_dietnum","water_dietnum",Poisson Regression
         # "alcohol_dietnum","yogurt_dietnum","dairy_dietnum","probiotic_dietnum","fruits_no_juice_dietnum",
         # "vegetables_dietnum","beans_soy_dietnum","whole_grains_dietnum","starch_dietnum"
         # "eggs_dietnum"
         # "processed_meat_dietnum","red_meat_dietnum","white_meat_dietnum","shellfish_dietnum","fish_dietnum",
         # "sweets_dietnum"
         )

  ### Specify a set of covariates
  # covar <- c("consent_age","sex",paste0("race",0:3),paste0("educ",0:7))
  covar <- c("consent_age")

  # Test when it works -----
  ## Test ZIP_qgcomp with covariates ----
  results <- ZIP_qgcomp(formatted_data = formatted_data,
                           x = x,
                           covar = covar,
                           q = 4)

  testthat::expect_equal(object = ncol(results), expected = 9)


  ## Test ZIP_qgcomp without covariates ----
  results <- ZIP_qgcomp(formatted_data = formatted_data,
                           covar = NULL,
                           x = x,
                           q = 4)

  testthat::expect_equal(object = ncol(results), expected = 9)


  ## Test ZIP_qgcomp with q = NULL ----
  results <- ZIP_qgcomp(formatted_data = formatted_data,
                         covar = covar,
                         x = x,
                         q = NULL)

  testthat::expect_equal(object = ncol(results), expected = 9)
})
