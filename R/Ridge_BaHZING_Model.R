#' Ridge_BaHZING_Model Function
#' This function implements the BaHZING model for microbiome data analysis.
#' @import R2jags
#' @import rjags
#' @import pscl
#' @import dplyr
#' @import tidyr
#' @import phyloseq
#' @import stringr
#' @import tibble
#' @importFrom utils globalVariables
#' @importFrom stats quantile update
#' @importFrom bayestestR p_direction p_rope p_map
#' @param formatted_data An object containing formatted microbiome data.
#' @param x A vector of column names of the exposures.
#' @param covar An optional vector of the column names of covariates.
#' @param n.chains An optional integer specifying the number of parallel chains
#' for the model in jags.model function. Default is 3.
#' @param n.adapt An optional integer specifying the number of iterations for
#' adaptation in jags.model function. Default is 100.
#' @param n.iter.burnin An optional integer specifying number of iterations in
#' update function. Default is 1,000.
#' @param n.iter.sample An optional integer specifying the number of iterations
#' in coda.samples function. Default is 5,000.
#' @param exposure_standardization Method for standardizing the exposures.
#' Should be one of "standard_normal" (the default), "quantile", or "none". If
#' "none", exposures are not standardized before analysis, and counterfactual
#' profiles must be specified by the user.
#' @param counterfactual_profiles A 2xP matrix or a vector with length of 2; P
#' is the number of exposures in x. If a 2xP matrix is provided,
#' the effect estimates for the mixture are interpreted as the estimated change
#' in the outcome when changing each exposure p in 1:P is changed from
#' `counterfactual_profiles[1,p]` to `counterfactual_profiles[2,p]`. If a vector of
#' length 2 is provided, the effect estimates for the mixture are interpreted as
#' the estimated change in the outcome when changing each exposure from
#' `counterfactual_profiles[1]` to `counterfactual_profiles[2]`. If
#' exposure_standardization = "standard_normal", then the default is c(-0.5, 0.5),
#' and the effect estimate is calculated based on increasing all exposures in
#' the mixture by one standard deviation. If exposure_standardization = "quantile",
#' then the default is c(0,1), and the effect estimate is calculated based on
#' increasing all exposures in the mixture by one quantile (where the number of
#' quantiles is based on the parameter q).
#' @param q An integer specifying the number of quantiles. Only required if
#' exposure_standardization = "quantile". If exposure_standardization =
#' "quantile" and q is not specified, then a default of q = 4 is used.
#' @param verbose If TRUE (default), function returns information a data quality
#' check.
#' @param return_all_estimates If FALSE (default), results do not include
#' the dispersion estimates from the BaHZING model.
#' @param ROPE_range Region of practical equivalence (ROPE) for calculating
#' p_rope. Default is c(-0.1, 0.1).
#' @return A data frame containing results of the Bayesian analysis, with the
#' following columns:
#' - taxa_full: Full Taxa information, including all levels of the taxonomy.
#' Taxanomic levels are split by two underscores ('__').
#' - taxa_name: Taxa name, which is the last level of the taxonomy.
#' - domain: domain of the taxa.
#' - exposure: Exposure name (either one of  the individual exposures, or the
#' mixture).
#' - component: Zero inflated model estimate or the Count model estimate.
#' - estimate: Point estimate of the posterior distributions.
#' - bci_lcl: 95% Bayesian Credible Interval Lower Limit. Calculated as the
#' equal tailed interval of posterior distributions using the quantiles method.
#' - bci_ucl: 95% Bayesian Credible Interval Upper Limit. Calculated as the
#' equal tailed interval of posterior distributions using the quantiles method.
#' - p_direction: The Probability of Direction, calculated with `bayestestR`. A
#' higher value suggests a higher probability that the estimate is strictly
#' positive or negative. In other words, the closer the value to 1, the higher
#' the probability that the estimate is non-zero. Values can not be less than
#' 50%. From `bayestestR`: also known as the Maximum Probability of Effect
#' (MPE). This can be interpreted as the probability that a parameter (described
#' by its posterior distribution) is strictly positive or negative (whichever
#' is the most probable). Although differently expressed, this index is fairly
#' similar (i.e., is strongly correlated) to the frequentist p-value.
#' - p_rope: The probability that the estimate is not within the Region of
#' practical equivalence (ROPE), calculated with `bayestestR`. The proportion
#' of the whole posterior distribution that doesn't lie within the `ROPE_range`.
#' - p_map: Bayesian equivalent of the p-value, calculated with `bayestestR`.
#'  From `bayestestR`:  p_map is related to the odds that a parameter (described
#'  by its posterior distribution) has against the null hypothesis (h0) using
#'   Mills' (2014, 2017) Objective Bayesian Hypothesis Testing framework. It
#'   corresponds to the density value at the null (e.g., 0) divided by the
#'   density at the Maximum A Posteriori (MAP).
#' @export
#' @name Ridge_BaHZING_Model

# Declare global variables
globalVariables(c("LibrarySize", "X2.5.", "X97.5.", "Mean",
                  "Exposure.Index", "taxa_index", "taxa_full",
                  "component", "estimate", "bci_lcl", "bci_ucl",
                  "domain", "taxa_name", "pdir","prope","pmap",
                  "name"))

Ridge_BaHZING_Model <- function(formatted_data,
                                x,
                                covar = NULL,
                                exposure_standardization = NULL,
                                n.chains = 3,
                                n.adapt = 100,
                                n.iter.burnin = 1000,
                                n.iter.sample = 5000,
                                counterfactual_profiles = NULL,
                                q = NULL,
                                verbose = TRUE,
                                return_all_estimates = FALSE,
                                ROPE_range = c(-0.1, 0.1)) {

  # 1. Check input data ----
  # Extract metadata file from formatted data
  exposure_covar_dat <- data.frame(formatted_data$Table)

  # Create covariate dataframe
  if (!is.null(covar)){
    W <- data.frame(exposure_covar_dat[covar])
    Q <- ncol(W)
  }

  # Create exposure dataframe
  if(!all(x %in% colnames(exposure_covar_dat))) {
    stop("Not all exposured are found in the formatted data")
  }
  X <- exposure_covar_dat[x]
  P <- ncol(X)

  # Set exposure_standardization if missing
  if(is.null(exposure_standardization)) {
    exposure_standardization = "standard_normal"
  }
  # Check exposure_standardization is within the bounds
  if(!(exposure_standardization %in% c("standard_normal", "quantile", "none"))){
    stop("exposure_standardization must be either standard_normal, quantile, or none")
  }

  # Give warning if exposure_standardization == "standard_normal" and q is specified
  if(exposure_standardization == "standard_normal" & !is.null(q)){
    message("Note: q is not required when exposure_standardization is standard_normal. q will be ignored.")
  }

  # If standardization is "none", then check to make sure that counterfactual_profiles is specified
  if(exposure_standardization  == "none" & is.null(counterfactual_profiles)){
    stop("counterfactual_profiles must be speficied if exposure_standardization is none.")
  }

  # 2. Counterfactual profiles -----------------------------------
  ## 2.a Set counterfactual_profiles if missing --------
  if(is.null(counterfactual_profiles)) {
    if(exposure_standardization == "standard_normal") {
      counterfactual_profiles <- c(-0.5,0.5)
    } else {
      counterfactual_profiles <- c(0,1)
    }
  }

  ## 2.b Check counterfactual profiles structure ----
  if (is.matrix(counterfactual_profiles)) { # Checks if matrix
    if (nrow(counterfactual_profiles) != 2) {
      stop("counterfactual_profiles must have 2 rows when provided as a matrix.")
    }
    if (ncol(counterfactual_profiles) != P) {
      stop(paste0("When provided as a matrix, the number of columns in counterfactual_profiles must be equal to the number of exposures in the model."))
    }
    if (!is.numeric(counterfactual_profiles)) {
      stop("counterfactual_profiles must be numeric.")
    }
  } else { # Checks if is vector
    if (is.vector(counterfactual_profiles)) {
      if (!is.numeric(counterfactual_profiles)) {
        stop("counterfactual_profiles must be numeric.")
      }
      if (length(counterfactual_profiles) != 2) {
        stop(paste0("counterfactual_profiles must have 2 elements when provided as a vector."))
      }
    } else {
      stop("counterfactual_profiles must be either a numeric 2xP matrix or a numeric vector with length P.")
    }
  }

  ## 2.c Set profiles --------------------------------
  if(is.matrix(counterfactual_profiles)) {
    profiles = counterfactual_profiles
  } else {
    profiles <- rbind(rep(counterfactual_profiles[1], P),
                      rep(counterfactual_profiles[2], P))
  }

  # Give warning if using quantiles but counterfactuals < 0
  if(exposure_standardization=="quantiles" &
     any(counterfactual_profiles<0 |
         counterfactual_profiles > q+1)){
    warning("Note: Quantiles are used, but counterfactual_profiles includes values outside of range. Estimates will be calculated based on the specified counterfactual_profiles, but results may not be interpretable.")
  }


  # 3. Scale exposures --------------------------------
  # Set default q if not provided
  if(is.null(q) & exposure_standardization == "quantile") {
    q = 4
  }

  # If using quantiles, quantize X
  if(exposure_standardization=="quantile") {
    probs <- seq(0, 1, length.out = q + 1)
    X.q <- apply(X, 2, function(v) {
      cut(v, breaks = c(-Inf, quantile(v, probs = probs, include.lowest = FALSE)), labels = FALSE)
    })
  }

  #If not quantized and not standardized, scale X
  if(exposure_standardization=="standard_normal") {
    X.q <- scale(X)
  }

  #If none, no scaling
  if(exposure_standardization == "none") {
    X.q <- X
  }

  # 4. Format microbiome matricies ----
  #Create outcome dataframe
  Y <- exposure_covar_dat[, grep("k__", names(exposure_covar_dat))]
  N <- nrow(Y)
  R <- ncol(Y)
  #Genus
  Z.s.g <- t(formatted_data$Species.Genus.Matrix)
  GenusData <- as.matrix(Y) %*% Z.s.g %>% as.data.frame()
  Genus.R <- ncol(GenusData)

  #Family
  Z.g.f <- t(formatted_data$Genus.Family.Matrix)
  FamilyData <- as.matrix(GenusData) %*% Z.g.f %>% as.data.frame()
  Family.R <- ncol(FamilyData)

  #Order
  Z.f.o <- t(formatted_data$Family.Order.Matrix)
  OrderData <- as.matrix(FamilyData) %*% Z.f.o %>% as.data.frame()
  Order.R <- ncol(OrderData)

  #Class
  Z.o.c <- t(formatted_data$Order.Class.Matrix)
  ClassData <- as.matrix(OrderData) %*% Z.o.c %>% as.data.frame()
  Class.R <- ncol(ClassData)

  #Phylum
  Z.c.p <- t(formatted_data$Class.Phylum.Matrix)
  PhylumData <- as.matrix(ClassData) %*% Z.c.p %>% as.data.frame()
  Phylum.R <- ncol(PhylumData)

  ## Create Library Size Offset
  L <- exposure_covar_dat[, grep("k__", names(exposure_covar_dat))]
  L <- L %>%
    mutate(LibrarySize=rowSums(across(everything())))
  L <- L %>%
    select(LibrarySize)

  # 5. Return "Sanity" Messages ----
  if(verbose == TRUE){
    message("#### Checking input data ####")
    message("Exposure and Covariate Data:")
    message(paste0("- Total sample size: ", N))
    message(paste0("- Number of exposures: ", P))

    message("Microbiome Data:")
    message(paste0("- Number of unique species in data: ",  R))
    message(paste0("- Number of unique genus in data: ",  Genus.R))
    message(paste0("- Number of unique family in data: ", Family.R))
    message(paste0("- Number of unique order in data: ",  Order.R))
    message(paste0("- Number of unique class in data: ",  Class.R))
    message(paste0("- Number of unique phylum in data: ", Phylum.R))

    message("#### Running BaHZING with the following parameters #### ")
    if (exposure_standardization == "standard_normal"){
      message("Exposure standardization: Standard Normal")
    }
    if (exposure_standardization == "none"){
      message("Exposure standardization: None")
    }
    if (exposure_standardization == "quantile"){
      message(paste0("Exposure standardization: Quantiles, with q = ", q))
    }
  }

  # 6. Model ----
  ## a. With Covariates----
  Ridge_BHRM.microbiome <-
    "model {
    for(r in 1:P.s) { # for each feature (species, genus, family, etc.)
      for(i in 1:N) { # loop through individuals

        # zero inflated negative binomial
        Y[i,r] ~ dnegbin(mu[i,r], disp[r])
        mu[i,r] <- disp[r]/(disp[r]+(1-zero[i,r])*lambda[i,r]) - 0.000001*zero[i,r]

        # means component
        log(lambda[i,r]) <- alpha[r] + inprod(beta[r,1:P.e], X.q[i,1:P.e]) + inprod(delta[r, 1:Q], W[i,1:Q]) + log(L[i,1])

        # zero inflation component
        zero[i,r] ~ dbern(pi[i,r])
        logit(pi[i,r]) <- alpha.zero[r] + inprod(beta.zero[r,1:P.e], X.q[i,1:P.e]) + inprod(delta.zero[r, 1:Q], W[i,1:Q]) + log(L[i,1])
      }
      # prior on dispersion parameter
      disp[r] ~ dunif(0,50)

      # prior on intercept
      alpha[r] ~ dnorm(0, 1.0E-02) # means component
      alpha.zero[r] ~ dnorm(0, 1.0E-02) # zero inflation component

      # prior on covariate effects
      for(q in 1:Q) {
        delta[r,q] ~ dnorm(0, 1.0E-02) # means component
        delta.zero[r,q] ~ dnorm(0, 1.0E-02) # zero inflation component
      }

      # prior on exposure effects
      for(p in 1:P.e) {
        beta[r,p] ~ dnorm(0, tau[r]) # means component
        beta.zero[r,p] ~ dnorm(0, tau.zero[r]) # zero inflation component
      }

      # prior on precision for exposure effects
      # means component
      tau[r] <- 1/(sigma[r]*sigma[r])
      sigma[r] ~ dunif(0,3)

      # zero inflation component
      tau.zero[r] <- 1/(sigma.zero[r]*sigma.zero[r])
      sigma.zero[r] ~ dunif(0,3)

      # g-estimation
      # means component
      eta.low[r] <- inprod(beta[r,1:P.e], profiles[1,1:P.e])
      eta.high[r] <- inprod(beta[r,1:P.e], profiles[2,1:P.e])
      psi[r] <- eta.high[r]-eta.low[r]

      # zero inflation component
      eta.low.zero[r] <- inprod(beta.zero[r,1:P.e], profiles[1,1:P.e])
      eta.high.zero[r] <- inprod(beta.zero[r,1:P.e], profiles[2,1:P.e])
      psi.zero[r] <- eta.high.zero[r]-eta.low.zero[r]
    }
  }"

  ## b. Without Covariates ----
  Ridge_BHRM_no_covariates.microbiome <-
    "model {
    for(r in 1:P.s) { # for each feature (species, genus, family, etc.)
      for(i in 1:N) { # loop through individuals
        # zero inflated negative binomial

        Y[i,r] ~ dnegbin(mu[i,r], disp[r])
        mu[i,r] <- disp[r]/(disp[r]+(1-zero[i,r])*lambda[i,r]) - 0.000001*zero[i,r]

        # means component
        log(lambda[i,r]) <- alpha[r] + inprod(beta[r,1:P.e], X.q[i,1:P.e]) + log(L[i,1])

        # zero inflation component
        zero[i,r] ~ dbern(pi[i,r])
        logit(pi[i,r]) <- alpha.zero[r] + inprod(beta.zero[r,1:P.e], X.q[i,1:P.e]) + log(L[i,1])
      }
      # prior on dispersion parameter
      disp[r] ~ dunif(0,50)

      # prior on intercept
      alpha[r] ~ dnorm(0, 1.0E-02) # means component
      alpha.zero[r] ~ dnorm(0, 1.0E-02) # zero inflation component

      # prior on exposure effects
      for(p in 1:P.e) {
        beta[r,p] ~ dnorm(0, tau[r]) # means component
        beta.zero[r,p] ~ dnorm(0, tau.zero[r]) # zero inflation component
      }

      # prior on precision for exposure effects
      # means component
      tau[r] <- 1/(sigma[r]*sigma[r])
      sigma[r] ~ dunif(0,3)

      # zero inflation component
      tau.zero[r] <- 1/(sigma.zero[r]*sigma.zero[r])
      sigma.zero[r] ~ dunif(0,3)

      # g-estimation
      # means component
      eta.low[r] <- inprod(beta[r,1:P.e], profiles[1,1:P.e])
      eta.high[r] <- inprod(beta[r,1:P.e], profiles[2,1:P.e])
      psi[r] <- eta.high[r]-eta.low[r]

      # zero inflation component
      eta.low.zero[r] <- inprod(beta.zero[r,1:P.e], profiles[1,1:P.e])
      eta.high.zero[r] <- inprod(beta.zero[r,1:P.e], profiles[2,1:P.e])
      psi.zero[r] <- eta.high.zero[r]-eta.low.zero[r]
    }
  }"

  # 7. Run Model ----
  ## a. With Covariates----
  if (!is.null(covar)) {
    ### Run JAGs Estimation for species----
    jdata <- list(N=N, Y=Y, P.s= R, X.q=X.q, P.e = P, W = W, Q = Q,
                  profiles=profiles,
                  L=L)

    var.s <- c("beta", "beta.zero","psi","disp")
    model.fit <- jags.model(file=textConnection( Ridge_BHRM.microbiome), data=jdata,
                            n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s,
                              n.iter=n.iter.sample, thin=1, progress.bar="text")

    ### Summary result for species----
    ### Calculate Mean, SD, and quantiles
    r <- summary(model.fit)
    results_species <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    ### Calculate p-values
    post_dist <-  as.data.frame(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]]))]

    ### p_direction
    pdir <- apply(post_dist, 2, function(x){
      p_direction(x = x, threshold = 0.05) %>%
        as.numeric()})

    ### p_rope
    prope <- apply(post_dist, 2, function(x){
      p_rope(x = x, rope = ROPE_range)$p_ROPE})

    ### p_map
    pmap <- apply(post_dist, 2, function(x){
      p_map(x = x) %>%
        as.numeric()})

    # Get dataframe of p-values
    p_value_df <- data.frame(name = names(pdir),
                             pdir = pdir,
                             prope = prope,
                             pmap = pmap)

    ## Create "component" variable
    results_species <- results_species %>%
      mutate(
        component=case_when(
          grepl("zero",rownames(.)) ~ "Zero-inflation model coefficients",
          grepl("beta",rownames(.)) ~ "Count model coefficients",
          grepl("psi",rownames(.))  ~ "Count model coefficients",
          grepl("disp",rownames(.)) ~ "Dispersion",
          TRUE ~ "Other"))%>%
      # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
      rename(estimate = Mean,
             bci_lcl = X2.5.,
             bci_ucl = X97.5.)%>%
      rownames_to_column("name")%>%
      # Combine results df with p-value df
      tidylog::left_join( p_value_df, by = "name")%>%
      mutate(name = ifelse(grepl("beta", name)|grepl("psi", name),
                           paste0("Species.", name),
                           name)) %>%
      column_to_rownames("name")


    ### Run JAGs Estimation for Genus----
    jdata <- list(N=N, Y=GenusData, P.s=Genus.R, X.q=X.q, P.e = P, W = W, Q = Q,
                  profiles=profiles,
                  L=L)

    var.s <- c("beta", "beta.zero","psi","disp")
    model.fit <- jags.model(file=textConnection(Ridge_BHRM.microbiome), data=jdata,
                            n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s,
                              n.iter=n.iter.sample, thin=1, progress.bar="text")

    #### Summary result for Genus----
    ### Calculate Mean, SD, and quantiles
    r <- summary(model.fit)
    results_Genus <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    ### Calculate p-values
    post_dist <-  as.data.frame(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]]))]

    ### p_direction
    pdir <- apply(post_dist, 2, function(x){
      p_direction(x = x, threshold = 0.05) %>%
        as.numeric()})

    ### p_rope
    prope <- apply(post_dist, 2, function(x){
      p_rope(x = x, rope = ROPE_range)$p_ROPE})

    ### p_map
    pmap <- apply(post_dist, 2, function(x){
      p_map(x = x) %>%
        as.numeric()})

    # Get dataframe of p-values
    p_value_df <- data.frame(name = names(pdir),
                             pdir = pdir,
                             prope = prope,
                             pmap = pmap)

    ## Create "component" variable
    results_Genus <- results_Genus %>%
      mutate(
        component=case_when(
          grepl("zero",rownames(.)) ~ "Zero-inflation model coefficients",
          grepl("beta",rownames(.)) ~ "Count model coefficients",
          grepl("psi",rownames(.))  ~ "Count model coefficients",
          grepl("disp",rownames(.)) ~ "Dispersion",

          TRUE ~ "Other"))%>%
      # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
      rename(estimate = Mean,
             bci_lcl = X2.5.,
             bci_ucl = X97.5.)%>%
      rownames_to_column("name")%>%
      # Combine results df with p-value df
      tidylog::left_join( p_value_df, by = "name")%>%
      mutate(name = ifelse(grepl("beta", name)|grepl("psi", name),
                           paste0("Genus.", name),
                           name)) %>%
      column_to_rownames("name")



    ### Run JAGs Estimation for Family----
    jdata <- list(N = N, Y= FamilyData, P.s=Family.R, X.q=X.q, P.e = P, W = W, Q = Q,
                  profiles=profiles,
                  L = L)

    var.s <- c("beta", "beta.zero","psi","disp")
    model.fit <- jags.model(file=textConnection(Ridge_BHRM.microbiome), data=jdata,
                            n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s,
                              n.iter=n.iter.sample, thin=1, progress.bar="text")

    #### Summary result for Family----
    ### Calculate Mean, SD, and quantiles
    r <- summary(model.fit)
    results_Family <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    ### Calculate p-values
    post_dist <-  as.data.frame(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]]))]

    ### p_direction
    pdir <- apply(post_dist, 2, function(x){
      p_direction(x = x, threshold = 0.05) %>%
        as.numeric()})

    ### p_rope
    prope <- apply(post_dist, 2, function(x){
      p_rope(x = x, rope = ROPE_range)$p_ROPE})

    ### p_map
    pmap <- apply(post_dist, 2, function(x){
      p_map(x = x) %>%
        as.numeric()})

    # Get dataframe of p-values
    p_value_df <- data.frame(name = names(pdir),
                             pdir = pdir,
                             prope = prope,
                             pmap = pmap)

    ## Create "component" variable
    results_Family <- results_Family %>%
      mutate(
        component=case_when(
          grepl("zero",rownames(.)) ~ "Zero-inflation model coefficients",
          grepl("beta",rownames(.)) ~ "Count model coefficients",
          grepl("psi",rownames(.))  ~ "Count model coefficients",
          grepl("disp",rownames(.)) ~ "Dispersion",

          TRUE ~ "Other"))%>%
      # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
      rename(estimate = Mean,
             bci_lcl = X2.5.,
             bci_ucl = X97.5.)%>%
      rownames_to_column("name")%>%
      # Combine results df with p-value df
      tidylog::left_join( p_value_df, by = "name")%>%
      mutate(name = ifelse(grepl("beta", name)|grepl("psi", name),
                           paste0("Family.", name),
                           name)) %>%
      column_to_rownames("name")


    ### Run JAGs Estimation for Order----
    jdata <- list(N=N, Y=OrderData, P.s=Order.R, X.q=X.q, P.e = P,W = W, Q = Q,
                  profiles=profiles,
                  L=L)

    var.s <- c("beta", "beta.zero","psi","disp")
    model.fit <- jags.model(file=textConnection( Ridge_BHRM.microbiome), data=jdata,
                            n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s,
                              n.iter=n.iter.sample, thin=1, progress.bar="text")

    #### Summary result for Order----
    ### Calculate Mean, SD, and quantiles
    r <- summary(model.fit)
    results_Order <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    ### Calculate p-values
    post_dist <-  as.data.frame(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]]))]

    ### p_direction
    pdir <- apply(post_dist, 2, function(x){
      p_direction(x = x, threshold = 0.05) %>%
        as.numeric()})

    ### p_rope
    prope <- apply(post_dist, 2, function(x){
      p_rope(x = x, rope = ROPE_range)$p_ROPE})

    ### p_map
    pmap <- apply(post_dist, 2, function(x){
      p_map(x = x) %>%
        as.numeric()})

    # Get dataframe of p-values
    p_value_df <- data.frame(name = names(pdir),
                             pdir = pdir,
                             prope = prope,
                             pmap = pmap)

    ## Create "component" variable
    results_Order <- results_Order %>%
      mutate(
        component=case_when(
          grepl("zero",rownames(.)) ~ "Zero-inflation model coefficients",
          grepl("beta",rownames(.)) ~ "Count model coefficients",
          grepl("psi",rownames(.))  ~ "Count model coefficients",
          grepl("disp",rownames(.)) ~ "Dispersion",

          TRUE ~ "Other"))%>%
      # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
      rename(estimate = Mean,
             bci_lcl = X2.5.,
             bci_ucl = X97.5.)%>%
      rownames_to_column("name")%>%
      # Combine results df with p-value df
      tidylog::left_join( p_value_df, by = "name")%>%
      mutate(name = ifelse(grepl("beta", name)|grepl("psi", name),
                           paste0("Order.", name),
                           name)) %>%
      column_to_rownames("name")


    ### Run JAGs Estimation for Class----
    jdata <- list(N=N, Y=ClassData, P.s=Class.R, X.q=X.q, P.e = P,W = W, Q = Q,
                  profiles=profiles,
                  L=L)

    var.s <- c("beta", "beta.zero","psi","disp")
    model.fit <- jags.model(file=textConnection( Ridge_BHRM.microbiome), data=jdata,
                            n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s,
                              n.iter=n.iter.sample, thin=1, progress.bar="text")

    #### Summary result for Class----
    ### Calculate Mean, SD, and quantiles
    r <- summary(model.fit)
    results_Class <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    ### Calculate p-values
    post_dist <-  as.data.frame(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]]))]

    ### p_direction
    pdir <- apply(post_dist, 2, function(x){
      p_direction(x = x, threshold = 0.05) %>%
        as.numeric()})

    ### p_rope
    prope <- apply(post_dist, 2, function(x){
      p_rope(x = x, rope = ROPE_range)$p_ROPE})

    ### p_map
    pmap <- apply(post_dist, 2, function(x){
      p_map(x = x) %>%
        as.numeric()})

    # Get dataframe of p-values
    p_value_df <- data.frame(name = names(pdir),
                             pdir = pdir,
                             prope = prope,
                             pmap = pmap)

    ## Create "component" variable
    results_Class <- results_Class %>%
      mutate(
        component=case_when(
          grepl("zero",rownames(.)) ~ "Zero-inflation model coefficients",
          grepl("beta",rownames(.)) ~ "Count model coefficients",
          grepl("psi",rownames(.))  ~ "Count model coefficients",
          grepl("disp",rownames(.)) ~ "Dispersion",

          TRUE ~ "Other"))%>%
      # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
      rename(estimate = Mean,
             bci_lcl = X2.5.,
             bci_ucl = X97.5.)%>%
      rownames_to_column("name")%>%
      # Combine results df with p-value df
      tidylog::left_join( p_value_df, by = "name")%>%
      mutate(name = ifelse(grepl("beta", name)|grepl("psi", name),
                           paste0("Class.", name),
                           name)) %>%
      column_to_rownames("name")



    ### Run JAGs Estimation for Phylum----
    jdata <- list(N=N, Y=PhylumData, P.s=Phylum.R, X.q=X.q, P.e = P, W = W, Q = Q,
                  profiles=profiles,
                  L=L)

    var.s <- c("beta", "beta.zero","psi","disp")
    model.fit <- jags.model(file=textConnection( Ridge_BHRM.microbiome), data=jdata,
                            n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s,
                              n.iter=n.iter.sample, thin=1, progress.bar="text")

    #### Summary result for Phylum----
    ### Calculate Mean, SD, and quantiles
    r <- summary(model.fit)
    results_Phylum <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    ### Calculate p-values
    post_dist <-  as.data.frame(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]]))]

    ### p_direction
    pdir <- apply(post_dist, 2, function(x){
      p_direction(x = x, threshold = 0.05) %>%
        as.numeric()})

    ### p_rope
    prope <- apply(post_dist, 2, function(x){
      p_rope(x = x, rope = ROPE_range)$p_ROPE})

    ### p_map
    pmap <- apply(post_dist, 2, function(x){
      p_map(x = x) %>%
        as.numeric()})

    # Get dataframe of p-values
    p_value_df <- data.frame(name = names(pdir),
                             pdir = pdir,
                             prope = prope,
                             pmap = pmap)

    ## Create "component" variable
    results_Phylum <- results_Phylum %>%
      mutate(
        component=case_when(
          grepl("zero",rownames(.)) ~ "Zero-inflation model coefficients",
          grepl("beta",rownames(.)) ~ "Count model coefficients",
          grepl("psi",rownames(.))  ~ "Count model coefficients",
          grepl("disp",rownames(.)) ~ "Dispersion",

          TRUE ~ "Other"))%>%
      # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
      rename(estimate = Mean,
             bci_lcl = X2.5.,
             bci_ucl = X97.5.)%>%
      rownames_to_column("name")%>%
      # Combine results df with p-value df
      tidylog::left_join( p_value_df, by = "name")%>%
      mutate(name = ifelse(grepl("beta", name)|grepl("psi", name),
                           paste0("Phylum.", name),
                           name)) %>%
      column_to_rownames("name")

    ### Combine result for Species, Genus, Family, Order, Class, Phylum-----
    results <- results_species %>%
      bind_rows(results_Genus)%>%
      bind_rows(results_Family) %>%
      bind_rows(results_Order) %>%
      bind_rows(results_Class) %>%
      bind_rows(results_Phylum)
  }else{
    ## b. Without covariates----
    ### Run JAGs Estimation for species----
    jdata <- list(N=N, Y=Y, P.s=R, X.q=X.q, P.e = P,
                  profiles=profiles,
                  L=L)

    var.s <- c("beta", "beta.zero","psi","disp")
    model.fit <- jags.model(file=textConnection( Ridge_BHRM_no_covariates.microbiome), data=jdata,
                            n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s,
                              n.iter=n.iter.sample, thin=1, progress.bar="text")

    ### Summary result for species----
    ### Calculate Mean, SD, and quantiles
    r <- summary(model.fit)
    results_species <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    ### Calculate p-values
    post_dist <-  as.data.frame(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]]))]

    ### p_direction
    pdir <- apply(post_dist, 2, function(x){
      p_direction(x = x, threshold = 0.05) %>%
        as.numeric()})

    ### p_rope
    prope <- apply(post_dist, 2, function(x){
      p_rope(x = x, rope = ROPE_range)$p_ROPE})

    ### p_map
    pmap <- apply(post_dist, 2, function(x){
      p_map(x = x) %>%
        as.numeric()})

    # Get dataframe of p-values
    p_value_df <- data.frame(name = names(pdir),
                             pdir = pdir,
                             prope = prope,
                             pmap = pmap)

    ## Create "component" variable
    results_species <- results_species %>%
      mutate(
        component=case_when(
          grepl("zero",rownames(.)) ~ "Zero-inflation model coefficients",
          grepl("beta",rownames(.)) ~ "Count model coefficients",
          grepl("psi",rownames(.))  ~ "Count model coefficients",
          grepl("disp",rownames(.)) ~ "Dispersion",

          TRUE ~ "Other"))%>%
      # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
      rename(estimate = Mean,
             bci_lcl = X2.5.,
             bci_ucl = X97.5.)%>%
      rownames_to_column("name")%>%
      # Combine results df with p-value df
      tidylog::left_join( p_value_df, by = "name")%>%
      mutate(name = ifelse(grepl("beta", name)|grepl("psi", name),
                           paste0("Species.", name),
                           name)) %>%
      column_to_rownames("name")


    ### Run JAGs Estimation for Genus----
    jdata <- list(N=N, Y=GenusData, P.s=Genus.R, X.q=X.q, P.e = P,
                  profiles=profiles,
                  L=L)

    var.s <- c("beta", "beta.zero","psi","disp")
    model.fit <- jags.model(file=textConnection( Ridge_BHRM_no_covariates.microbiome), data=jdata,
                            n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s,
                              n.iter=n.iter.sample, thin=1, progress.bar="text")

    #### Summary result for Genus----
    ### Calculate Mean, SD, and quantiles
    r <- summary(model.fit)
    results_Genus <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    ### Calculate p-values
    post_dist <-  as.data.frame(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]]))]

    ### p_direction
    pdir <- apply(post_dist, 2, function(x){
      p_direction(x = x, threshold = 0.05) %>%
        as.numeric()})

    ### p_rope
    prope <- apply(post_dist, 2, function(x){
      p_rope(x = x, rope = ROPE_range)$p_ROPE})

    ### p_map
    pmap <- apply(post_dist, 2, function(x){
      p_map(x = x) %>%
        as.numeric()})

    # Get dataframe of p-values
    p_value_df <- data.frame(name = names(pdir),
                             pdir = pdir,
                             prope = prope,
                             pmap = pmap)

    ## Create "component" variable
    results_Genus <- results_Genus %>%
      mutate(
        component=case_when(
          grepl("zero",rownames(.)) ~ "Zero-inflation model coefficients",
          grepl("beta",rownames(.)) ~ "Count model coefficients",
          grepl("psi",rownames(.))  ~ "Count model coefficients",
          grepl("disp",rownames(.)) ~ "Dispersion",

          TRUE ~ "Other"))%>%
      # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
      rename(estimate = Mean,
             bci_lcl = X2.5.,
             bci_ucl = X97.5.)%>%
      rownames_to_column("name")%>%
      # Combine results df with p-value df
      tidylog::left_join( p_value_df, by = "name")%>%
      mutate(name = ifelse(grepl("beta", name)|grepl("psi", name),
                           paste0("Genus.", name),
                           name)) %>%
      column_to_rownames("name")



    ### Run JAGs Estimation for Family----
    jdata <- list(N=N, Y=FamilyData, P.s=Family.R, X.q=X.q, P.e = P,
                  profiles=profiles,
                  L=L)

    var.s <- c("beta", "beta.zero","psi","disp")
    model.fit <- jags.model(file=textConnection( Ridge_BHRM_no_covariates.microbiome), data=jdata,
                            n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s,
                              n.iter=n.iter.sample, thin=1, progress.bar="text")

    #### Summary result for Family----
    ### Calculate Mean, SD, and quantiles
    r <- summary(model.fit)
    results_Family <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    ### Calculate p-values
    post_dist <-  as.data.frame(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]]))]

    ### p_direction
    pdir <- apply(post_dist, 2, function(x){
      p_direction(x = x, threshold = 0.05) %>%
        as.numeric()})

    ### p_rope
    prope <- apply(post_dist, 2, function(x){
      p_rope(x = x, rope = ROPE_range)$p_ROPE})

    ### p_map
    pmap <- apply(post_dist, 2, function(x){
      p_map(x = x) %>%
        as.numeric()})

    # Get dataframe of p-values
    p_value_df <- data.frame(name = names(pdir),
                             pdir = pdir,
                             prope = prope,
                             pmap = pmap)

    ## Create "component" variable
    results_Family <- results_Family %>%
      mutate(
        component=case_when(
          grepl("zero",rownames(.)) ~ "Zero-inflation model coefficients",
          grepl("beta",rownames(.)) ~ "Count model coefficients",
          grepl("psi",rownames(.))  ~ "Count model coefficients",
          grepl("disp",rownames(.)) ~ "Dispersion",

          TRUE ~ "Other"))%>%
      # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
      rename(estimate = Mean,
             bci_lcl = X2.5.,
             bci_ucl = X97.5.)%>%
      rownames_to_column("name")%>%
      # Combine results df with p-value df
      tidylog::left_join( p_value_df, by = "name")%>%
      mutate(name = ifelse(grepl("beta", name)|grepl("psi", name),
                           paste0("Family.", name),
                           name)) %>%
      column_to_rownames("name")


    ### Run JAGs Estimation for Order----
    jdata <- list(N=N, Y=OrderData, P.s=Order.R, X.q=X.q, P.e = P,
                  profiles=profiles,
                  L=L)

    var.s <- c("beta", "beta.zero","psi","disp")
    model.fit <- jags.model(file=textConnection( Ridge_BHRM_no_covariates.microbiome), data=jdata,
                            n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s,
                              n.iter=n.iter.sample, thin=1, progress.bar="text")

    #### Summary result for Order----
    ### Calculate Mean, SD, and quantiles
    r <- summary(model.fit)
    results_Order <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    ### Calculate p-values
    post_dist <-  as.data.frame(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]]))]

    ### p_direction
    pdir <- apply(post_dist, 2, function(x){
      p_direction(x = x, threshold = 0.05) %>%
        as.numeric()})

    ### p_rope
    prope <- apply(post_dist, 2, function(x){
      p_rope(x = x, rope = ROPE_range)$p_ROPE})

    ### p_map
    pmap <- apply(post_dist, 2, function(x){
      p_map(x = x) %>%
        as.numeric()})

    # Get dataframe of p-values
    p_value_df <- data.frame(name = names(pdir),
                             pdir = pdir,
                             prope = prope,
                             pmap = pmap)

    ## Create "component" variable
    results_Order <- results_Order %>%
      mutate(
        component=case_when(
          grepl("zero",rownames(.)) ~ "Zero-inflation model coefficients",
          grepl("beta",rownames(.)) ~ "Count model coefficients",
          grepl("psi",rownames(.))  ~ "Count model coefficients",
          grepl("disp",rownames(.)) ~ "Dispersion",

          TRUE ~ "Other"))%>%
      # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
      rename(estimate = Mean,
             bci_lcl = X2.5.,
             bci_ucl = X97.5.)%>%
      rownames_to_column("name")%>%
      # Combine results df with p-value df
      tidylog::left_join( p_value_df, by = "name")%>%
      mutate(name = ifelse(grepl("beta", name)|grepl("psi", name),
                           paste0("Order.", name),
                           name)) %>%
      column_to_rownames("name")


    ### Run JAGs Estimation for Class----
    jdata <- list(N=N, Y=ClassData, P.s=Class.R, X.q=X.q, P.e = P,
                  profiles=profiles,
                  L=L)

    var.s <- c("beta", "beta.zero","psi","disp")
    model.fit <- jags.model(file=textConnection( Ridge_BHRM_no_covariates.microbiome), data=jdata,
                            n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s,
                              n.iter=n.iter.sample, thin=1, progress.bar="text")

    #### Summary result for Class----
    ### Calculate Mean, SD, and quantiles
    r <- summary(model.fit)
    results_Class <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    ### Calculate p-values
    post_dist <-  as.data.frame(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]]))]

    ### p_direction
    pdir <- apply(post_dist, 2, function(x){
      p_direction(x = x, threshold = 0.05) %>%
        as.numeric()})

    ### p_rope
    prope <- apply(post_dist, 2, function(x){
      p_rope(x = x, rope = ROPE_range)$p_ROPE})

    ### p_map
    pmap <- apply(post_dist, 2, function(x){
      p_map(x = x) %>%
        as.numeric()})

    # Get dataframe of p-values
    p_value_df <- data.frame(name = names(pdir),
                             pdir = pdir,
                             prope = prope,
                             pmap = pmap)

    ## Create "component" variable
    results_Class <- results_Class %>%
      mutate(
        component=case_when(
          grepl("zero",rownames(.)) ~ "Zero-inflation model coefficients",
          grepl("beta",rownames(.)) ~ "Count model coefficients",
          grepl("psi",rownames(.))  ~ "Count model coefficients",
          grepl("disp",rownames(.)) ~ "Dispersion",

          TRUE ~ "Other"))%>%
      # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
      rename(estimate = Mean,
             bci_lcl = X2.5.,
             bci_ucl = X97.5.)%>%
      rownames_to_column("name")%>%
      # Combine results df with p-value df
      tidylog::left_join( p_value_df, by = "name")%>%
      mutate(name = ifelse(grepl("beta", name)|grepl("psi", name),
                           paste0("Class.", name),
                           name)) %>%
      column_to_rownames("name")



    ### Run JAGs Estimation for Phylum----
    jdata <- list(N=N, Y=PhylumData, P.s=Phylum.R, X.q=X.q, P.e = P,
                  profiles=profiles,
                  L=L)

    var.s <- c("beta", "beta.zero","psi","disp")
    model.fit <- jags.model(file=textConnection( Ridge_BHRM_no_covariates.microbiome), data=jdata,
                            n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s,
                              n.iter=n.iter.sample, thin=1, progress.bar="text")

    #### Summary result for Phylum----
    ### Calculate Mean, SD, and quantiles
    r <- summary(model.fit)
    results_Phylum <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    ### Calculate p-values
    post_dist <-  as.data.frame(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]]))]

    ### p_direction
    pdir <- apply(post_dist, 2, function(x){
      p_direction(x = x, threshold = 0.05) %>%
        as.numeric()})

    ### p_rope
    prope <- apply(post_dist, 2, function(x){
      p_rope(x = x, rope = ROPE_range)$p_ROPE})

    ### p_map
    pmap <- apply(post_dist, 2, function(x){
      p_map(x = x) %>%
        as.numeric()})

    # Get dataframe of p-values
    p_value_df <- data.frame(name = names(pdir),
                             pdir = pdir,
                             prope = prope,
                             pmap = pmap)

    ## Create "component" variable
    results_Phylum <- results_Phylum %>%
      mutate(
        component=case_when(
          grepl("zero",rownames(.)) ~ "Zero-inflation model coefficients",
          grepl("beta",rownames(.)) ~ "Count model coefficients",
          grepl("psi",rownames(.))  ~ "Count model coefficients",
          grepl("disp",rownames(.)) ~ "Dispersion",

          TRUE ~ "Other"))%>%
      # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
      rename(estimate = Mean,
             bci_lcl = X2.5.,
             bci_ucl = X97.5.)%>%
      rownames_to_column("name")%>%
      # Combine results df with p-value df
      tidylog::left_join( p_value_df, by = "name")%>%
      mutate(name = ifelse(grepl("beta", name)|grepl("psi", name),
                           paste0("Phylum.", name),
                           name)) %>%
      column_to_rownames("name")

    ### Combine result for Species, Genus, Family, Order, Class, Phylum-----
    results <- results_species %>%
      bind_rows(results_Genus)%>%
      bind_rows(results_Family) %>%
      bind_rows(results_Order) %>%
      bind_rows(results_Class) %>%
      bind_rows(results_Phylum)
  }

  # 8. Format output data----
  phylum   <- colnames(PhylumData)
  class    <- colnames(ClassData)
  order    <- colnames(OrderData)
  family   <- colnames(FamilyData)
  genus    <- colnames(GenusData)
  species  <- colnames(Y)
  exposure <- colnames(X)

  results <- results %>%
    mutate(taxa_index = str_remove(rownames(results),"..*\\[") %>%
             str_remove(.,",.*$") %>%
             str_remove(.,"].*$") %>% as.numeric(),
           Exposure.Index = str_remove(rownames(results),"..*\\,"))%>%
    mutate(Exposure.Index=ifelse(
      grepl("disp",Exposure.Index)|
        grepl("psi",Exposure.Index),
      NA,Exposure.Index)) %>%
    mutate(Exposure.Index = suppressWarnings(str_remove(Exposure.Index,"]") %>%
                                               as.numeric())) %>%
    mutate(
      taxa_full=case_when(
        grepl("Phylum", rownames(results)) ~ paste0(phylum[taxa_index]),
        grepl("Class"  ,rownames(results)) ~ paste0(class[taxa_index]),
        grepl("Order"  ,rownames(results)) ~ paste0(order[taxa_index]),
        grepl("Family" ,rownames(results)) ~ paste0(family[taxa_index]),
        grepl("Genus"  ,rownames(results)) ~ paste0(genus[taxa_index]),
        grepl("Species",rownames(results)) ~ paste0(species[taxa_index])),
      domain=case_when(
        grepl("Phylum" ,rownames(results))  ~ "Phylum",
        grepl("Class"  ,rownames(results))   ~ "Class",
        grepl("Order"  ,rownames(results))   ~ "Order",
        grepl("Family" ,rownames(results))  ~ "Family",
        grepl("Genus"  ,rownames(results))   ~ "Genus",
        grepl("Species",rownames(results)) ~ "Species"),
      exposure=paste0(exposure[Exposure.Index])) %>%
    # Modify Exposure variable
    mutate(exposure=case_when(
      grepl("psi",rownames(.)) ~ "Mixture",
      grepl("disp",rownames(.)) ~ "Dispersion",

      TRUE ~ exposure))

  # Get Taxa and domain information
  results2 <- results %>%
    mutate(taxa_full= ifelse(grepl("disp", rownames(results)),paste0(species[taxa_index]),taxa_full),
           taxa_name = sub(".*__", "", taxa_full),
           domain = ifelse(exposure == "Dispersion",
                           "Species", domain))

  # Remove "disp" and "omega" estimates
  if(!return_all_estimates){
    results2 <- results2 %>%
      filter(!grepl("disp",rownames(results2)))
  }

  # Remove rownames
  rownames(results2) <- NULL

  # Select final variables
  results2 <- results2 %>%
    select(taxa_full, taxa_name, domain, exposure,component,
           estimate,bci_lcl,bci_ucl,pdir,prope,pmap)

  return(results2)
  }

