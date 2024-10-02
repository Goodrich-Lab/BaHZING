test_that("test Format_BaHZING", {
  PS <- iHMP
  formatted_data <- Format_BaHZING(PS)
  # Test format data
  testthat::expect_equal(object = length(formatted_data), expected = 6)
  testthat::expect_true(is.list(formatted_data))
  testthat::expect_true("Table" %in% names(formatted_data))
})

# test_that("taxa_are_rows transposes correctly",{
#   # OTU table with taxa are in rows
#   PS <- iHMP
#   otu_matrix <- iHMP@otu_table@.Data %>% t()
#   OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)
#   PS@otu_table <- OTU
#
#   formatted_data <- Format_BaHZING(PS)
#   # Check that taxa are in rows
#   testthat::expect_equal(object = length(formatted_data), expected = 6)
#   testthat::expect_true(is.list(formatted_data))
#   testthat::expect_true("Table" %in% names(formatted_data))
# })

test_that("Error thrown when taxonomic levels are less than 2", {

  PS <- iHMP
  # Modify the taxonomic table to have only one level
  tax_table(PS) <- tax_table(tax_table(PS)[,1])

  # Expect that the function stops with an error message
  testthat::expect_error(Format_BaHZING(PS), "Need > 1 taxonomic level")
  })

# test_that("Error thrown when missing Kingdom-level information", {
#   PS <- iHMP
#
#   # Set one Kingdom level to NA
#   tax_table(PS)[1, "Kingdom"] <- NA
#
#   # Expect an error to be thrown due to missing Kingdom-level information
#   testthat::expect_error(Format_BaHZING(PS),
#                          "Missing Kindom-level information for an ASV/OTU. Remove unidentified bacteria and rerun.")
#
# })

test_that("If species level not present, create species column", {
  PS <- iHMP
  # phyloseq object without a 'Species' column
  tax_table(PS) <- tax_table(PS)[ , !colnames(tax_table(PS)) %in% "Species"]

  formatted_data <- Format_BaHZING(PS)

  # Extract the taxonomic table from the result
  taxa_table_result <- formatted_data$Table[[1]]

  # Check if 'Species' column is created
  testthat::expect_true(grepl("s__unclassified",
                              colnames(taxa_table_result)[ncol(taxa_table_result)]))
  testthat::expect_true(grepl("s__unclassified",
                              colnames(taxa_table_result)[ncol(taxa_table_result)-10]))
  testthat::expect_true(grepl("s__unclassified",
                              colnames(taxa_table_result)[ncol(taxa_table_result)-20]))
})

