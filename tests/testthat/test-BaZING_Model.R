test_that("test BaZING_Model", {
  PS <- simulated_data
  formatted_data <- Format_BaZING(PS)
  covar <- c("energy","age","sex","race0","race1","race2","education0","education1","education2")
  x <- c("pfda","pfna","pfhxs","pfhps","pfpes","pfos","pfoa")
  results <- BaZING_Model(formatted_data, covar, x)
  # Test format data
  testthat::expect_equal(object = ncol(results), expected = 7)
})


