#' BaZING_Model Function
#' This function implements the BaZING model for microbiome data analysis.
#'
#' Load Packages
#' @import R2jags
#' @import rjags
#' @import pscl
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import phyloseq
#' @import stringr

#' @param formatted_data An object containing formatted microbiome data.
#' @param covar A vector of covariates.
#' @param x A vector of exposures.
#' @param n.chain An integer specifying the number of parallel chains for the model in jags.model function
#' @param n.adapt An integer specifying the number of iterations for adaptation in jags.model function
#' @param n.iter.update An integer specifying number of iterations in update function
#' @param n.iter.coda An integer specifying the number of iterations in coda.samples function
#' @param standardized Logical. If TRUE, uses standardized exposures
#' @param counterfactual_profiles A 2xP matrix or a vector with length of 2; P is the column number of exposure dataframe.
#' @param q An integer specifying the number of quantile groups
#' @return A data frame containing results of the Bayesian analysis.


BaZING_Model <- function(formatted_data,
                         covar,
                         x,
                         exposure_standardization,
                         n.chains = 3,
                         n.adapt = 100,
                         n.iter.update = 2,
                         n.iter.coda = 2,
                         # standardized = FALSE,
                         counterfactual_profiles = c(-0.5, 0.5),
                         q = 4) {


  #Set default counterfactual profiles
  default <- c(-0.5,0.5)
  #Extract metaddata file from formatted data
  Object <- data.frame(formatted_data$Table)

  #Create covariate dataframe
  W <- data.frame(Object[covar])
  # W <- data.frame(Sex,Education)
  Q <- ncol(W)

  #Create exposure dataframe
  X <- Object[x]
  P <- ncol(X)

  if(is.matrix(counterfactual_profiles)){
    if(nrow(counterfactual_profiles)!=2|ncol(counterfactual_profiles)!=P){
      "Counterfactural Profiles should be a 2xP matrix or a vector with length 2"
    }else{
      profiles = counterfactual_profiles
    }
  }else{
    if(length(counterfactual_profiles) != 2){
      "Counterfactural Profiles should be a 2xP matrix or a vector with length 2"
    }
    if(exposure_standardization=="quantiles" & all(counterfactual_profiles==default)){
      counterfactual_profiles = c(0, 1)
    }
    profiles <- rbind(rep(counterfactual_profiles[1], P), rep(counterfactual_profiles[2], P))
  }

  if(!(exposure_standardization %in% c("standard_normal", "quantile"))){
    stop("exposure_standardization must be either standard_normal or quantile")
  }

  #If quantiles specified as true, quantize X
  if (exposure_standardization=="quantile") {
    probs <- seq(0, 1, length.out = q + 1)
    X.q <- apply(X, 2, function(v) {
      cut(v, breaks = quantile(v, probs = probs, include.lowest = TRUE), labels = FALSE)
    })
  }

  #If not quantized and not standardized, scale X
  if(exposure_standardization=="standard_normal") {
    X.q <- scale(X)
  # } else {
  #   X.q <- X
  }

  # Create profiles matrix
  #Create outcome dataframe
  Y <- Object[, grep("k__", names(Object))]
  N <- nrow(Y)
  R <- ncol(Y)

  #Format microbiome matricies
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

  #Create Library Size Offset
  L <- Object[, grep("k__", names(Object))]
  L <- L %>%
    mutate(LibrarySize=rowSums(across(everything())))
  L <- L %>%
    select(LibrarySize)

  # Hierarchical Model ----
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
  # model.fit <- jags.model(file=textConnection(BHRM.microbiome), data=jdata, n.chains=1, n.adapt=100, quiet=F)
  # update(model.fit, n.iter=10, progress.bar="text")
  # model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=50, thin=1, progress.bar="text")
  model.fit <- jags.model(file=textConnection(BHRM.microbiome), data=jdata, n.chains=n.chains, n.adapt=n.adapt, quiet=F)
  update(model.fit, n.iter=n.iter.update, progress.bar="text")
  model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=n.iter.coda, thin=1, progress.bar="text")
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
