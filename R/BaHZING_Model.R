#' BaHZING_Model Function
#' This function implements the BaHZING model for microbiome data analysis.
#'
#' Load Packages
#' @import R2jags
#' @import rjags
#' @import pscl
#' @import dplyr
#' @import tidyr
#' @import phyloseq
#' @import stringr
#' @importFrom stats quantile update
#' @param formatted_data An object containing formatted microbiome data.
#' @param x A vector of column names of the exposures.
#' @param covar An optional vector of the column names of covariates.
#' @param n.chains An optional integer specifying the number of parallel chains
#' for the model in jags.model function. Default is 3.
#' @param n.adapt An optional integer specifying the number of iterations for
#' adaptation in jags.model function. Default is 5000.
#' @param n.iter.burnin An optional integer specifying number of iterations in
#' update function. Default is 10000.
#' @param n.iter.sample An optional integer specifying the number of iterations
#' in coda.samples function. Default is 10000.
#' @param exposure_standardization Method for standardizing the exposures.
#' Should be one of "standard_normal" (the default), "quantile", or "none". If
#' "none", exposures are not standardized before analysis, and counterfactual
#' profiles must be specified by the user.
#' @param counterfactual_profiles A 2xP matrix or a vector with length of 2; P
#' is the number of exposures in x. If a 2xP matrix is provided,
#' the effect estimates for the mixture are interpreted as the estimated change
#' in the outcome when changing each exposure p in 1:P is changed from
#' counterfactual_profiles[1,p] to counterfactual_profiles[2,p]. If a vector of
#' length 2 is provided, the effect estimates for the mixture are interpreted as
#' the estimated change in the outcome when changing each exposure from
#' counterfactual_profiles[1] to counterfactual_profiles[2]. If
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
#' @return A data frame containing results of the Bayesian analysis.

#' @export
#' @name BaHZING_Model

# Declare global variables
utils::globalVariables(c("LibrarySize", "X2.5.", "X97.5.", "Mean",
                         "Exposure.Index", "Taxa.Index", "Taxa",
                         "Component", "OR", "OR.ll", "OR.ul", "sig", "Domain"))

BaHZING_Model <- function(formatted_data,
                          x,
                          covar = NULL,
                          exposure_standardization = NULL,
                          n.chains = 3,
                          n.adapt = 5000,
                          n.iter.burnin = 10000,
                          n.iter.sample = 10000,
                          counterfactual_profiles = NULL,
                          q = NULL,
                          verbose = TRUE) {

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

  # Give warning if using qualtiles but counterfactuals < 0
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

  ## Create Library Size Offset
  L <- exposure_covar_dat[, grep("k__", names(exposure_covar_dat))]
  L <- L %>%
    mutate(LibrarySize=rowSums(across(everything())))
  L <- L %>%
    select(LibrarySize)

  # 3. Return "Sanity" Messages ----
  if(verbose == TRUE){
    message("#### Checking input data ####")
    message("Exposure and Covariate Data:")
    message(paste0("- Total sample size: ", N))
    message(paste0("- Number of exposures: ", P))

    message("Microbiome Data:")
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

  # 4. Run Model ----
  if (!is.null(covar)) {
    # Hierarchical Model with covariates----
    BHRM.microbiome <-
      "model {
    for(r in 1:R) {
      for(i in 1:N) {
        Y[i,r] ~ dnegbin(mu[i,r], disp[r])
        mu[i,r] <- disp[r]/(disp[r]+(1-zero[i,r])*lambda[i,r]) - 0.000001*zero[i,r]
        log(lambda[i,r]) <- alpha[r] + inprod(species.beta[r,1:P], X.q[i,1:P]) + inprod(delta[r, 1:Q], W[i,1:Q]) + log(L[i,1])

        # zero-inflation
        zero[i,r] ~ dbern(pi[i,r])
        logit(pi[i,r]) <- alpha.zero[r] + inprod(species.beta.zero[r,1:P], X.q[i,1:P]) + inprod(delta.zero[r, 1:Q], W[i,1:Q]) + log(L[i,1])
      }
      # prior on dispersion parameter
      disp[r] ~ dunif(0,50)

      # prior on intercept
      alpha[r] ~ dnorm(0, 1.0E-02)
      alpha.zero[r] ~ dnorm(0, 1.0E-02)

      # prior on proportion of non-zeros
      omega[r] ~ dunif(0,1)

      # prior on covariate effects
      for(q in 1:Q) {
        delta[r,q] ~ dnorm(0, 1.0E-02)
        delta.zero[r,q] ~ dnorm(0, 1.0E-02)
      }

      # prior on exposure effects
      for(p in 1:P) {
        # species.beta.zero[r,p] ~ dnorm(0, 1.0E-02)
        species.beta[r,p] ~ dnorm(mu.species[r,p], tau[r])
        mu.species[r,p] <- inprod(genus.beta[1:Genus.R,p], GenusData[r,1:Genus.R])
        #Zero inflation component
        species.beta.zero[r,p] ~ dnorm(mu.species.zero[r,p], tau[r])
        mu.species.zero[r,p] <- inprod(genus.beta.zero[1:Genus.R,p], GenusData[r,1:Genus.R])
      }

      # prior on precision
      tau[r] <- 1/(sigma[r]*sigma[r])
      sigma[r] ~ dunif(0,3)

      # g-estimation
      species.eta.low[r] <- inprod(species.beta[r,1:P], profiles[1,1:P])
      species.eta.high[r] <- inprod(species.beta[r,1:P], profiles[2,1:P])
      species.psi[r] <- species.eta.high[r]-species.eta.low[r]
      # zero-inflation
      species.eta.low.zero[r] <- inprod(species.beta.zero[r,1:P], profiles[1,1:P])
      species.eta.high.zero[r] <- inprod(species.beta.zero[r,1:P], profiles[2,1:P])
      species.psi.zero[r] <- species.eta.high.zero[r]-species.eta.low.zero[r]
    }

    # Genus level
    for(g.r in 1:Genus.R) {
      for(p in 1:P) {
        genus.beta[g.r,p] ~ dnorm(mu.family[g.r,p],genus.tau[g.r]) # should this be shared effects across exposures at genus level? or shared effects across all genus by exposure?
        mu.family[g.r,p] <- inprod(family.beta[1:Family.R,p], FamilyData[g.r,1:Family.R])
        #Zero inflation component
        genus.beta.zero[g.r,p] ~ dnorm(mu.family.zero[g.r,p],genus.tau[g.r]) # should this be shared effects across exposures at genus level? or shared effects across all genus by exposure?
        mu.family.zero[g.r,p] <- inprod(family.beta.zero[1:Family.R,p], FamilyData[g.r,1:Family.R])
      }
      # prior on precision
      genus.tau[g.r] <- 1/(genus.sigma[g.r]*genus.sigma[g.r])
      genus.sigma[g.r] ~ dunif(0,3)

      # g-estimation
      genus.eta.low[g.r] <- inprod(genus.beta[g.r,1:P], profiles[1,1:P])
      genus.eta.high[g.r] <- inprod(genus.beta[g.r,1:P], profiles[2,1:P])
      genus.psi[g.r] <- genus.eta.high[g.r]-genus.eta.low[g.r]
      #zero inflation
      genus.eta.low.zero[g.r] <- inprod(genus.beta.zero[g.r,1:P], profiles[1,1:P])
      genus.eta.high.zero[g.r] <- inprod(genus.beta.zero[g.r,1:P], profiles[2,1:P])
      genus.psi.zero[g.r] <- genus.eta.high.zero[g.r]-genus.eta.low.zero[g.r]
    }

    # Family level
    for(f.r in 1:Family.R) {
      for(p in 1:P) {
        family.beta[f.r,p] ~ dnorm(mu.order[f.r,p], family.tau[f.r])
        mu.order[f.r,p] <- inprod(order.beta[1:Order.R,p], OrderData[f.r,1:Order.R])
        #Zero inflation component
        family.beta.zero[f.r,p] ~ dnorm(mu.order.zero[f.r,p], family.tau[f.r])
        mu.order.zero[f.r,p] <- inprod(order.beta.zero[1:Order.R,p], OrderData[f.r,1:Order.R])

      }
      # prior on precision
      family.tau[f.r] <- 1/(family.sigma[f.r]*family.sigma[f.r])
      family.sigma[f.r] ~ dunif(0,3)

      # g-estimation
      family.eta.low[f.r] <- inprod(family.beta[f.r,1:P], profiles[1,1:P])
      family.eta.high[f.r] <- inprod(family.beta[f.r,1:P], profiles[2,1:P])
      family.psi[f.r] <- family.eta.high[f.r]-family.eta.low[f.r]
      #zero inflation
      family.eta.low.zero[f.r] <- inprod(family.beta.zero[f.r,1:P], profiles[1,1:P])
      family.eta.high.zero[f.r] <- inprod(family.beta.zero[f.r,1:P], profiles[2,1:P])
      family.psi.zero[f.r] <- family.eta.high.zero[f.r]-family.eta.low.zero[f.r]
    }

    # Order level
    for(o.r in 1:Order.R) {
      for(p in 1:P) {
        order.beta[o.r,p] ~ dnorm(mu.class[o.r,p], order.tau[o.r])
        mu.class[o.r,p] <- inprod(class.beta[1:Class.R,p], ClassData[o.r,1:Class.R])
        #Zero inflation component
        order.beta.zero[o.r,p] ~ dnorm(mu.class.zero[o.r,p], order.tau[o.r])
        mu.class.zero[o.r,p] <- inprod(class.beta.zero[1:Class.R,p], ClassData[o.r,1:Class.R])
      }
      # prior on precision
      order.tau[o.r] <- 1/(order.sigma[o.r]*order.sigma[o.r])
      order.sigma[o.r] ~ dunif(0,3)

      # g-estimation
      order.eta.low[o.r] <- inprod(order.beta[o.r,1:P], profiles[1,1:P])
      order.eta.high[o.r] <- inprod(order.beta[o.r,1:P], profiles[2,1:P])
      order.psi[o.r] <- order.eta.high[o.r]-order.eta.low[o.r]
      #zero infl
      order.eta.low.zero[o.r] <- inprod(order.beta.zero[o.r,1:P], profiles[1,1:P])
      order.eta.high.zero[o.r] <- inprod(order.beta.zero[o.r,1:P], profiles[2,1:P])
      order.psi.zero[o.r] <- order.eta.high.zero[o.r]-order.eta.low.zero[o.r]
    }

    # Class level
    for(c.r in 1:Class.R) {
      for(p in 1:P) {
        class.beta[c.r,p] ~ dnorm(mu.phylum[c.r,p], class.tau[c.r])
        mu.phylum[c.r,p] <- inprod(phylum.beta[1:Phylum.R,p], PhylumData[c.r,1:Phylum.R])
        #Zero inflation component
        class.beta.zero[c.r,p] ~ dnorm(mu.phylum.zero[c.r,p], class.tau[c.r])
        mu.phylum.zero[c.r,p] <- inprod(phylum.beta.zero[1:Phylum.R,p], PhylumData[c.r,1:Phylum.R])
      }
      # prior on precision
      class.tau[c.r] <- 1/(class.sigma[c.r]*class.sigma[c.r])
      class.sigma[c.r] ~ dunif(0,3)

      # g-estimation
      class.eta.low[c.r] <- inprod(class.beta[c.r,1:P], profiles[1,1:P])
      class.eta.high[c.r] <- inprod(class.beta[c.r,1:P], profiles[2,1:P])
      class.psi[c.r] <- class.eta.high[c.r]-class.eta.low[c.r]

      #zero component
      class.eta.low.zero[c.r] <- inprod(class.beta.zero[c.r,1:P], profiles[1,1:P])
      class.eta.high.zero[c.r] <- inprod(class.beta.zero[c.r,1:P], profiles[2,1:P])
      class.psi.zero[c.r] <- class.eta.high.zero[c.r]-class.eta.low.zero[c.r]
    }

    # Phylum level
    for(p.r in 1:Phylum.R) {
      for(p in 1:P) {
        phylum.beta[p.r,p] ~ dnorm(0, phylum.tau[p.r])
        #Zero inflation component
        phylum.beta.zero[p.r,p] ~ dnorm(0, phylum.tau[p.r])
      }
      # prior on precision
      phylum.tau[p.r] <- 1/(phylum.sigma[p.r]*phylum.sigma[p.r])
      phylum.sigma[p.r] ~ dunif(0,3)

      # g-estimation
      phylum.eta.low[p.r] <- inprod(phylum.beta[p.r,1:P], profiles[1,1:P])
      phylum.eta.high[p.r] <- inprod(phylum.beta[p.r,1:P], profiles[2,1:P])
      phylum.psi[p.r] <- phylum.eta.high[p.r]-phylum.eta.low[p.r]

      #Zero inflation
      phylum.eta.low.zero[p.r] <- inprod(phylum.beta.zero[p.r,1:P], profiles[1,1:P])
      phylum.eta.high.zero[p.r] <- inprod(phylum.beta.zero[p.r,1:P], profiles[2,1:P])
      phylum.psi.zero[p.r] <- phylum.eta.high.zero[p.r]-phylum.eta.low.zero[p.r]
    }

  }"

    # Run JAGs Estimation
    jdata <- list(N=N, Y=Y, R=R, X.q=X.q, W=W, P=P, Q=Q,
                  GenusData=GenusData, Genus.R=Genus.R,
                  Family.R=Family.R, FamilyData=FamilyData,
                  Order.R=Order.R, OrderData=OrderData,
                  Class.R=Class.R, ClassData=ClassData,
                  Phylum.R=Phylum.R, PhylumData=PhylumData,
                  profiles=profiles,L=L)
    var.s <- c("species.beta", "genus.beta", "family.beta", "order.beta", "class.beta", "phylum.beta", "species.beta.zero", "genus.beta.zero", "family.beta.zero", "order.beta.zero", "class.beta.zero", "phylum.beta.zero","species.psi","genus.psi","family.psi","order.psi","class.psi","phylum.psi","species.psi.zero","genus.psi.zero","family.psi.zero","order.psi.zero","class.psi.zero","phylum.psi.zero", "omega","disp")
    model.fit <- jags.model(file=textConnection(BHRM.microbiome), data=jdata, n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=n.iter.sample, thin=1, progress.bar="text")
    # summarize results
    r <- summary(model.fit)
    results <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    #Calculate significance based on Bayesian Interval
    results <- results %>%
      mutate(sig = ifelse((X2.5.<0 & X97.5.<0) | (X2.5.>0 & X97.5.>0), "*","N.S."))
    results <- results %>%
      mutate(OR=exp(Mean),
             OR.ll=exp(X2.5.),
             OR.ul=exp(X97.5.))
    #Format output names
    phylum <- colnames(PhylumData)
    class <- colnames(ClassData)
    order <- colnames(OrderData)
    family <- colnames(FamilyData)
    genus <- colnames(GenusData)
    species <- colnames(Y)
    Exposure <- colnames(X)
    results$Taxa.Index <- str_remove(rownames(results),"..*\\[")
    results$Taxa.Index <- str_remove(results$Taxa.Index,",.*$")
    results$Taxa.Index <- str_remove(results$Taxa.Index,"]")
    results$Taxa.Index <- as.numeric(results$Taxa.Index)
    results$Exposure.Index <- str_remove(rownames(results),"..*\\,")
    results <- results %>%
      mutate(Exposure.Index=ifelse(grepl("disp",Exposure.Index)|grepl("omega",Exposure.Index)|grepl("psi",Exposure.Index),NA,Exposure.Index))
    results$Exposure.Index <- str_remove(results$Exposure.Index,"]")
    results$Exposure.Index <- as.numeric(results$Exposure.Index)

    results <- results %>%
      mutate(Taxa=case_when(grepl("phylum",rownames(results)) ~ paste0(phylum[Taxa.Index]),
                            grepl("class",rownames(results)) ~ paste0(class[Taxa.Index]),
                            grepl("order",rownames(results)) ~ paste0(order[Taxa.Index]),
                            grepl("family",rownames(results)) ~ paste0(family[Taxa.Index]),
                            grepl("genus",rownames(results)) ~ paste0(genus[Taxa.Index]),
                            grepl("species",rownames(results)) ~ paste0(species[Taxa.Index])),
             Exposure=paste0(Exposure[Exposure.Index]))
    results <- results %>%
      mutate(Component=ifelse(grepl("zero",rownames(results)),"Means","Probability"))

    results <- results %>%
      mutate(Exposure=ifelse(grepl("disp",rownames(results)),paste0("Dispersion"),
                             ifelse(grepl("omega",rownames(results)),paste0("Omega"),Exposure)))
    results <- results %>%
      mutate(Exposure=ifelse(grepl("psi",rownames(results)),"Mixture",Exposure))
    results <- results %>%
      mutate(Taxa=ifelse(grepl("disp",rownames(results)),paste0(species[Taxa.Index]),Taxa),
             Taxa=ifelse(grepl("omega",rownames(results)),paste0(species[Taxa.Index]),Taxa))
    results <- results %>%
      select(Taxa,Exposure,Component,OR,OR.ll,OR.ul,sig)

    return(results)
  } else {
    # Hierarchical Model without Covariates----
    BHRM.microbiome <-
      "model {
    for(r in 1:R) {
      for(i in 1:N) {
        Y[i,r] ~ dnegbin(mu[i,r], disp[r])
        mu[i,r] <- disp[r]/(disp[r]+(1-zero[i,r])*lambda[i,r]) - 0.000001*zero[i,r]
        log(lambda[i,r]) <- alpha[r] + inprod(species.beta[r,1:P], X.q[i,1:P]) + log(L[i,1])

        # zero-inflation
        zero[i,r] ~ dbern(pi[i,r])
        logit(pi[i,r]) <- alpha.zero[r] + inprod(species.beta.zero[r,1:P], X.q[i,1:P]) + log(L[i,1])
      }
      # prior on dispersion parameter
      disp[r] ~ dunif(0,50)

      # prior on intercept
      alpha[r] ~ dnorm(0, 1.0E-02)
      alpha.zero[r] ~ dnorm(0, 1.0E-02)

      # prior on proportion of non-zeros
      omega[r] ~ dunif(0,1)

      # # prior on covariate effects
      # for(q in 1:Q) {
      #   delta[r,q] ~ dnorm(0, 1.0E-02)
      #   delta.zero[r,q] ~ dnorm(0, 1.0E-02)
      # }

      # prior on exposure effects
      for(p in 1:P) {
        # species.beta.zero[r,p] ~ dnorm(0, 1.0E-02)
        species.beta[r,p] ~ dnorm(mu.species[r,p], tau[r])
        mu.species[r,p] <- inprod(genus.beta[1:Genus.R,p], GenusData[r,1:Genus.R])
        #Zero inflation component
        species.beta.zero[r,p] ~ dnorm(mu.species.zero[r,p], tau[r])
        mu.species.zero[r,p] <- inprod(genus.beta.zero[1:Genus.R,p], GenusData[r,1:Genus.R])
      }

      # prior on precision
      tau[r] <- 1/(sigma[r]*sigma[r])
      sigma[r] ~ dunif(0,3)

      # g-estimation
      species.eta.low[r] <- inprod(species.beta[r,1:P], profiles[1,1:P])
      species.eta.high[r] <- inprod(species.beta[r,1:P], profiles[2,1:P])
      species.psi[r] <- species.eta.high[r]-species.eta.low[r]
      # zero-inflation
      species.eta.low.zero[r] <- inprod(species.beta.zero[r,1:P], profiles[1,1:P])
      species.eta.high.zero[r] <- inprod(species.beta.zero[r,1:P], profiles[2,1:P])
      species.psi.zero[r] <- species.eta.high.zero[r]-species.eta.low.zero[r]
    }

    # Genus level
    for(g.r in 1:Genus.R) {
      for(p in 1:P) {
        genus.beta[g.r,p] ~ dnorm(mu.family[g.r,p],genus.tau[g.r]) # should this be shared effects across exposures at genus level? or shared effects across all genus by exposure?
        mu.family[g.r,p] <- inprod(family.beta[1:Family.R,p], FamilyData[g.r,1:Family.R])
        #Zero inflation component
        genus.beta.zero[g.r,p] ~ dnorm(mu.family.zero[g.r,p],genus.tau[g.r]) # should this be shared effects across exposures at genus level? or shared effects across all genus by exposure?
        mu.family.zero[g.r,p] <- inprod(family.beta.zero[1:Family.R,p], FamilyData[g.r,1:Family.R])
      }
      # prior on precision
      genus.tau[g.r] <- 1/(genus.sigma[g.r]*genus.sigma[g.r])
      genus.sigma[g.r] ~ dunif(0,3)

      # g-estimation
      genus.eta.low[g.r] <- inprod(genus.beta[g.r,1:P], profiles[1,1:P])
      genus.eta.high[g.r] <- inprod(genus.beta[g.r,1:P], profiles[2,1:P])
      genus.psi[g.r] <- genus.eta.high[g.r]-genus.eta.low[g.r]
      #zero inflation
      genus.eta.low.zero[g.r] <- inprod(genus.beta.zero[g.r,1:P], profiles[1,1:P])
      genus.eta.high.zero[g.r] <- inprod(genus.beta.zero[g.r,1:P], profiles[2,1:P])
      genus.psi.zero[g.r] <- genus.eta.high.zero[g.r]-genus.eta.low.zero[g.r]
    }

    # Family level
    for(f.r in 1:Family.R) {
      for(p in 1:P) {
        family.beta[f.r,p] ~ dnorm(mu.order[f.r,p], family.tau[f.r])
        mu.order[f.r,p] <- inprod(order.beta[1:Order.R,p], OrderData[f.r,1:Order.R])
        #Zero inflation component
        family.beta.zero[f.r,p] ~ dnorm(mu.order.zero[f.r,p], family.tau[f.r])
        mu.order.zero[f.r,p] <- inprod(order.beta.zero[1:Order.R,p], OrderData[f.r,1:Order.R])

      }
      # prior on precision
      family.tau[f.r] <- 1/(family.sigma[f.r]*family.sigma[f.r])
      family.sigma[f.r] ~ dunif(0,3)

      # g-estimation
      family.eta.low[f.r] <- inprod(family.beta[f.r,1:P], profiles[1,1:P])
      family.eta.high[f.r] <- inprod(family.beta[f.r,1:P], profiles[2,1:P])
      family.psi[f.r] <- family.eta.high[f.r]-family.eta.low[f.r]
      #zero inflation
      family.eta.low.zero[f.r] <- inprod(family.beta.zero[f.r,1:P], profiles[1,1:P])
      family.eta.high.zero[f.r] <- inprod(family.beta.zero[f.r,1:P], profiles[2,1:P])
      family.psi.zero[f.r] <- family.eta.high.zero[f.r]-family.eta.low.zero[f.r]
    }

    # Order level
    for(o.r in 1:Order.R) {
      for(p in 1:P) {
        order.beta[o.r,p] ~ dnorm(mu.class[o.r,p], order.tau[o.r])
        mu.class[o.r,p] <- inprod(class.beta[1:Class.R,p], ClassData[o.r,1:Class.R])
        #Zero inflation component
        order.beta.zero[o.r,p] ~ dnorm(mu.class.zero[o.r,p], order.tau[o.r])
        mu.class.zero[o.r,p] <- inprod(class.beta.zero[1:Class.R,p], ClassData[o.r,1:Class.R])
      }
      # prior on precision
      order.tau[o.r] <- 1/(order.sigma[o.r]*order.sigma[o.r])
      order.sigma[o.r] ~ dunif(0,3)

      # g-estimation
      order.eta.low[o.r] <- inprod(order.beta[o.r,1:P], profiles[1,1:P])
      order.eta.high[o.r] <- inprod(order.beta[o.r,1:P], profiles[2,1:P])
      order.psi[o.r] <- order.eta.high[o.r]-order.eta.low[o.r]
      #zero infl
      order.eta.low.zero[o.r] <- inprod(order.beta.zero[o.r,1:P], profiles[1,1:P])
      order.eta.high.zero[o.r] <- inprod(order.beta.zero[o.r,1:P], profiles[2,1:P])
      order.psi.zero[o.r] <- order.eta.high.zero[o.r]-order.eta.low.zero[o.r]
    }

    # Class level
    for(c.r in 1:Class.R) {
      for(p in 1:P) {
        class.beta[c.r,p] ~ dnorm(mu.phylum[c.r,p], class.tau[c.r])
        mu.phylum[c.r,p] <- inprod(phylum.beta[1:Phylum.R,p], PhylumData[c.r,1:Phylum.R])
        #Zero inflation component
        class.beta.zero[c.r,p] ~ dnorm(mu.phylum.zero[c.r,p], class.tau[c.r])
        mu.phylum.zero[c.r,p] <- inprod(phylum.beta.zero[1:Phylum.R,p], PhylumData[c.r,1:Phylum.R])
      }
      # prior on precision
      class.tau[c.r] <- 1/(class.sigma[c.r]*class.sigma[c.r])
      class.sigma[c.r] ~ dunif(0,3)

      # g-estimation
      class.eta.low[c.r] <- inprod(class.beta[c.r,1:P], profiles[1,1:P])
      class.eta.high[c.r] <- inprod(class.beta[c.r,1:P], profiles[2,1:P])
      class.psi[c.r] <- class.eta.high[c.r]-class.eta.low[c.r]

      #zero component
      class.eta.low.zero[c.r] <- inprod(class.beta.zero[c.r,1:P], profiles[1,1:P])
      class.eta.high.zero[c.r] <- inprod(class.beta.zero[c.r,1:P], profiles[2,1:P])
      class.psi.zero[c.r] <- class.eta.high.zero[c.r]-class.eta.low.zero[c.r]
    }

    # Phylum level
    for(p.r in 1:Phylum.R) {
      for(p in 1:P) {
        phylum.beta[p.r,p] ~ dnorm(0, phylum.tau[p.r])
        #Zero inflation component
        phylum.beta.zero[p.r,p] ~ dnorm(0, phylum.tau[p.r])
      }
      # prior on precision
      phylum.tau[p.r] <- 1/(phylum.sigma[p.r]*phylum.sigma[p.r])
      phylum.sigma[p.r] ~ dunif(0,3)

      # g-estimation
      phylum.eta.low[p.r] <- inprod(phylum.beta[p.r,1:P], profiles[1,1:P])
      phylum.eta.high[p.r] <- inprod(phylum.beta[p.r,1:P], profiles[2,1:P])
      phylum.psi[p.r] <- phylum.eta.high[p.r]-phylum.eta.low[p.r]

      #Zero inflation
      phylum.eta.low.zero[p.r] <- inprod(phylum.beta.zero[p.r,1:P], profiles[1,1:P])
      phylum.eta.high.zero[p.r] <- inprod(phylum.beta.zero[p.r,1:P], profiles[2,1:P])
      phylum.psi.zero[p.r] <- phylum.eta.high.zero[p.r]-phylum.eta.low.zero[p.r]
    }

  }"

    # Run JAGs Estimation
    jdata <- list(N=N, Y=Y, R=R, X.q=X.q, P=P,
                  GenusData=GenusData, Genus.R=Genus.R,
                  Family.R=Family.R, FamilyData=FamilyData,
                  Order.R=Order.R, OrderData=OrderData,
                  Class.R=Class.R, ClassData=ClassData,
                  Phylum.R=Phylum.R, PhylumData=PhylumData,
                  profiles=profiles,L=L)
    var.s <- c("species.beta", "genus.beta", "family.beta", "order.beta",
               "class.beta", "phylum.beta", "species.beta.zero",
               "genus.beta.zero", "family.beta.zero", "order.beta.zero",
               "class.beta.zero", "phylum.beta.zero","species.psi","genus.psi",
               "family.psi","order.psi","class.psi","phylum.psi",
               "species.psi.zero","genus.psi.zero","family.psi.zero",
               "order.psi.zero","class.psi.zero","phylum.psi.zero",
               "omega","disp")
    model.fit <- jags.model(file=textConnection(BHRM.microbiome), data=jdata, n.chains=n.chains, n.adapt=n.adapt, quiet=F)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=n.iter.sample, thin=1, progress.bar="text")
    # summarize results
    r <- summary(model.fit)
    results <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

    #Calculate significance based on Bayesian Interval
    results <- results %>%
      mutate(sig = ifelse((X2.5.<0 & X97.5.<0) | (X2.5.>0 & X97.5.>0), "*","N.S."))
    results <- results %>%
      mutate(OR=exp(Mean),
             OR.ll=exp(X2.5.),
             OR.ul=exp(X97.5.))
    #Format output names
    phylum <- colnames(PhylumData)
    class <- colnames(ClassData)
    order <- colnames(OrderData)
    family <- colnames(FamilyData)
    genus <- colnames(GenusData)
    species <- colnames(Y)
    Exposure <- colnames(X)
    results$Taxa.Index <- str_remove(rownames(results),"..*\\[")
    results$Taxa.Index <- str_remove(results$Taxa.Index,",.*$")
    results$Taxa.Index <- str_remove(results$Taxa.Index,"]")
    results$Taxa.Index <- as.numeric(results$Taxa.Index)
    results$Exposure.Index <- str_remove(rownames(results),"..*\\,")
    results <- results %>%
      mutate(Exposure.Index=ifelse(grepl("disp",Exposure.Index)|grepl("omega",Exposure.Index)|grepl("psi",Exposure.Index),NA,Exposure.Index))
    results$Exposure.Index <- str_remove(results$Exposure.Index,"]")
    results$Exposure.Index <- as.numeric(results$Exposure.Index)

    results <- results %>%
      mutate(Taxa=case_when(grepl("phylum",rownames(results)) ~ paste0(phylum[Taxa.Index]),
                            grepl("class",rownames(results)) ~ paste0(class[Taxa.Index]),
                            grepl("order",rownames(results)) ~ paste0(order[Taxa.Index]),
                            grepl("family",rownames(results)) ~ paste0(family[Taxa.Index]),
                            grepl("genus",rownames(results)) ~ paste0(genus[Taxa.Index]),
                            grepl("species",rownames(results)) ~ paste0(species[Taxa.Index])),
             Exposure=paste0(Exposure[Exposure.Index]))

    results <- results %>%
      mutate(Component=ifelse(grepl("zero",rownames(results)),"Means","Probability"))

    results <- results %>%
      mutate(Exposure=ifelse(grepl("disp",rownames(results)),paste0("Dispersion"),
                             ifelse(grepl("omega",rownames(results)),paste0("Omega"),Exposure)))
    results <- results %>%
      mutate(Exposure=ifelse(grepl("psi",rownames(results)),"Mixture",Exposure))

    results <- results %>%
      mutate(Taxa=ifelse(grepl("disp",rownames(results)),paste0(species[Taxa.Index]),Taxa),
             Taxa=ifelse(grepl("omega",rownames(results)),paste0(species[Taxa.Index]),Taxa))
    results <- results %>%
      select(Taxa,Exposure,Component,OR,OR.ll,OR.ul,sig)

    return(results)
  }
}
