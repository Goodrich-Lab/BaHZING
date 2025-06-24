#' ZIP_qgcomp Function
#' This function implements the BaHZING model for microbiome data analysis.
#' @import R2jags
#' @import rjags
#' @import pscl
#' @import dplyr
#' @import tidyr
#' @import phyloseq
#' @import stringr
#' @import qgcomp
#' @importFrom utils globalVariables
#' @importFrom stats quantile update
#' @importFrom purrr map_dfr
#' @param formatted_data An object containing formatted microbiome data.
#' @param x A vector of column names of the exposures.
#' @param covar An optional vector of the column names of covariates.
#' @param q An integer specifying the number of quantiles.
#' @param verbose If TRUE (default), function returns information a data quality
#' check.
#' @return A data frame containing results of the qgcomp analysis, with the
#' following columns:
#' - taxa_full: Full Taxa information, including all levels of the taxonomy.
#' Taxanomic levels are split by two underscores ('__').
#' - taxa_name: Taxa name, which is the last level of the taxonomy.
#' - domain: domain of the taxa.
#' - exposure: Exposure name (either one of  the individual exposures, or the
#' mixture).
#' - estimate: Point estimate or mixture estimate.
#' - lcl: 95% Interval Lower Limit. Calculated as estimate - 1.96 × standard error
#' - ucl: 95% Interval Upper Limit. Calculated as estimate + 1.96 × standard error.
#' - p_value: P-value for the hypothesis that the effect estimate equals zero.
#' - model: Indicates the method used in the qgcomp() analysis.
#'   - "Poisson": Poisson regression model (from qgcomp()).
#' @export
#' @name ZIP_qgcomp
#'

ZIP_qgcomp <- function(formatted_data,
                        x,
                        covar = NULL,
                        q = NULL,
                        verbose = TRUE) {
  # 1. Check input data ----
  # Extract metadata file from formatted data
  exposure_covar_dat <- data.frame(formatted_data$Table)

  # Create covariate dataframe
  if (!is.null(covar)){
    W <- data.frame(exposure_covar_dat[covar])
    Q <- ncol(W)
  }else{
    Q = 0
  }

  # Create exposure dataframe
  if(!all(x %in% colnames(exposure_covar_dat))) {
    stop("Not all exposured are found in the formatted data")
  }
  X <- exposure_covar_dat[x]
  P <- ncol(X)


  #Create outcome dataframe
  Y.s <- exposure_covar_dat[, grep("k__", names(exposure_covar_dat))]

  #Genus
  GenusData <- as.data.frame(t(formatted_data$Species.Genus.Matrix))
  Genus.R <- ncol(GenusData)
  numGenusPerSpecies <- as.numeric(apply(GenusData, 1, sum))
  #Family
  FamilyData <- as.data.frame(t(formatted_data$Genus.Family.Matrix))
  Family.R <- ncol(FamilyData)
  numFamilyPerGenus <- as.numeric(apply(FamilyData, 1, sum))
  #Order
  OrderData <- as.data.frame(t(formatted_data$Family.Order.Matrix))
  Order.R <- ncol(OrderData)
  numOrderPerorder <- as.numeric(apply(OrderData, 1, sum))
  #Class
  ClassData <- as.data.frame(t(formatted_data$Order.Class.Matrix))
  Class.R <- ncol(ClassData)
  numClassPerOrder <- as.numeric(apply(ClassData, 1, sum))
  #Phylum
  PhylumData <- as.data.frame(t(formatted_data$Class.Phylum.Matrix))
  Phylum.R <- ncol(PhylumData)
  numPhylumPerClass <- as.numeric(apply(PhylumData, 1, sum))


  Y.g <- as.data.frame(as.matrix(Y.s)%*%as.matrix(GenusData))
  Y.f <- as.data.frame(as.matrix(Y.g)%*%as.matrix(FamilyData))
  Y.o <- as.data.frame(as.matrix(Y.f)%*%as.matrix(OrderData))
  Y.c <- as.data.frame(as.matrix(Y.o)%*%as.matrix(ClassData))
  Y.p <- as.data.frame(as.matrix(Y.c)%*%as.matrix(PhylumData))

  exposure_covar_dat <-  exposure_covar_dat %>%
    bind_cols(Y.g) %>%
    bind_cols(Y.f) %>%
    bind_cols(Y.o) %>%
    bind_cols(Y.c) %>%
    bind_cols(Y.p)
  # outcome including different domains
  Y <- exposure_covar_dat[, grep("k__", names(exposure_covar_dat))]
  N <- nrow(Y)
  R <- ncol(Y)

  # 2. Return "Sanity" Messages ----
  if(verbose == TRUE){
    message("#### Checking input data ####")
    message("Exposure and Covariate Data:")
    message(paste0("- Total sample size: ", N))
    message(paste0("- Number of exposures: ", P))
    message(paste0("- Number of covariates ", Q))

    message("Microbiome Data:")
    message(paste0("- Number of unique species in data: ",  dim(Y.s)[2]))
    message(paste0("- Number of unique genus in data: ",  Genus.R))
    message(paste0("- Number of unique family in data: ", Family.R))
    message(paste0("- Number of unique order in data: ",  Order.R))
    message(paste0("- Number of unique class in data: ",  Class.R))
    message(paste0("- Number of unique phylum in data: ", Phylum.R))

    message("#### Running qgcomp with the following parameters #### ")
    if (is.null(q)) {
      message("No quantile transformation applied.")
    } else {
      message(paste0("Quantiles, with q = ", q))
    }
  }

  ZING_Model <- function(r, dat, exposure, covariates, q ) {
      model_results <- data.frame()

      # Build the model formula
      formula_str <- paste0("`", r, "`", " ~ ", paste(c(exposure, covariates), collapse = " + "))
      m0 <- as.formula(formula_str)
      # Use Poisson
      mod1 <- qgcomp(f = m0, expnms = exposure, data = dat, q = q, family = poisson())
      coef_mat <- summary(mod1$fit)$coefficients
      mixture <- data.frame(
        taxa_full = r, exposure = "Mixture",
        estimate = coef_mat[2, 1],
        sd = coef_mat[2, 2],
        p_value = coef_mat[2, 4]
      )%>%
        mutate(
          lcl = estimate - 1.96*sd,
          ucl = estimate + 1.96*sd
        )

      individual <- data.frame(
        taxa_full = rep(r, length(exposure)),
        exposure = exposure,
        estimate = coef_mat[2:(length(exposure)+1), 1],
        sd = coef_mat[2:(length(exposure)+1), 2],
        p_value = coef_mat[2:(length(exposure)+1), 4]
      ) %>%
        mutate(
          lcl = estimate - 1.96*sd,
          ucl = estimate + 1.96*sd
        )

      model_results <- rbind(mixture, individual)
      model_results$model <- "Poisson"
      rownames(model_results) <- NULL

      return(model_results)
      }

    res <- map_dfr(colnames(Y),
                   ~ZING_Model(.x, dat = exposure_covar_dat, exposure = x, covariates = covar, q = q))

    # Adding domain
    res2 <- res %>% mutate(domain=case_when(grepl("_s__", taxa_full) ~ "Species",
                                            grepl("_g__", taxa_full)   ~ "Genus",
                                            grepl("_f__", taxa_full)  ~ "Family",
                                            grepl("_o__", taxa_full)   ~ "Order",
                                            grepl("_c__", taxa_full)   ~ "Class",
                                            grepl("_p__", taxa_full)  ~ "Phylum"),
                           taxa_name = sub(".*__", "", taxa_full)) %>%
      select(taxa_full, taxa_name, domain, exposure,
             estimate,lcl,ucl,p_value, model)

  return(res2)
}

