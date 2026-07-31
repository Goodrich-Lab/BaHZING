#' BaHZING_Model Function
#' This function implements the BaHZING model for microbiome data analysis.
#' @import R2jags
#' @import rjags
#' @import pscl
#' @import dplyr
#' @import tidyr
#' @import phyloseq
#' @import stringr
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
#' the dispersion and omega estimates from the BaHZING model.
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
#' @name BaHZING_Model

# Declare global variables
globalVariables(c("LibrarySize", "X2.5.", "X97.5.", "Mean",
                  "Exposure.Index", "taxa_index", "taxa_full",
                  "component", "estimate", "bci_lcl", "bci_ucl",
                  "domain", "taxa_name", "pdir","prope","pmap",
                  "name"))

# JAGS model text generators -------------------------------------------------
# The six taxonomic levels (species -> genus -> family -> order -> class ->
# phylum) share one recurring pattern: each level's exposure effect is drawn
# from a normal distribution centered on its parent level's corresponding
# effect (Equation 6), with species (the base/likelihood level) and phylum
# (the top, with no parent to borrow from) as the two special cases. These
# generators replace what used to be ~180 lines of hand-maintained, near-
# identical literal JAGS text per model variant (with/without covariates).
#
# All five generators below use glue::glue(.open="<<", .close=">>", .trim=FALSE):
# custom delimiters because JAGS's own "{ }" block syntax would otherwise
# collide with glue's default interpolation markers, and .trim=FALSE because
# glue's default (.trim=TRUE) strips leading whitespace and the trailing
# newline from its result - harmless for a block that ends the whole
# template, but it silently glues adjacent fragments onto the same line
# wherever one block's output is spliced into another (found the hard way:
# a comment ended up concatenated onto the same line as the code after it).

#' Build the precision-prior and g-estimation lines for one taxonomic level
#'
#' Generates the JAGS lines shared byte-for-byte by all six taxonomic levels:
#' the precision priors for both the count-model and zero-inflation dispersion
#' (`<level>.tau`/`<level>.sigma`, count and zero-inflation), followed by the
#' g-estimation block that computes each level's mixture contrast
#' (`<level>.psi`) from the low/high counterfactual `profiles` (Equation 7).
#' Used by every level via [.bahzing_level_block()] and directly by
#' [.bahzing_species_block()].
#'
#' @param level Character. The level's variable-name prefix, e.g. `"genus"`,
#'   `"species"`. Used both for the `<level>.tau`/`<level>.sigma` precision
#'   nodes and the `<level>.psi`/`<level>.eta.*` g-estimation nodes.
#' @param idx Character. The JAGS loop index variable for this level, e.g.
#'   `"g.r"` for genus, `"r"` for species.
#' @return A length-1 character string: the JAGS text for this block,
#'   ending in a trailing newline so it composes cleanly with adjacent
#'   fragments.
#' @keywords internal
#' @noRd
.bahzing_precision_gestimation_block <- function(level, idx) {
  tau_prefix <- paste0(level, ".")
  glue::glue(
"      # prior on precision
      <<tau_prefix>>tau[<<idx>>] <- 1/(<<tau_prefix>>sigma[<<idx>>]*<<tau_prefix>>sigma[<<idx>>])
      <<tau_prefix>>sigma[<<idx>>] ~ dunif(0,3)
      <<tau_prefix>>tau.zero[<<idx>>] <- 1/(<<tau_prefix>>sigma.zero[<<idx>>]*<<tau_prefix>>sigma.zero[<<idx>>])
      <<tau_prefix>>sigma.zero[<<idx>>] ~ dunif(0,3)

      # g-estimation
      <<level>>.eta.low[<<idx>>] <- inprod(<<level>>.beta[<<idx>>,1:P], profiles[1,1:P])
      <<level>>.eta.high[<<idx>>] <- inprod(<<level>>.beta[<<idx>>,1:P], profiles[2,1:P])
      <<level>>.psi[<<idx>>] <- <<level>>.eta.high[<<idx>>]-<<level>>.eta.low[<<idx>>]
      #zero inflation
      <<level>>.eta.low.zero[<<idx>>] <- inprod(<<level>>.beta.zero[<<idx>>,1:P], profiles[1,1:P])
      <<level>>.eta.high.zero[<<idx>>] <- inprod(<<level>>.beta.zero[<<idx>>,1:P], profiles[2,1:P])
      <<level>>.psi.zero[<<idx>>] <- <<level>>.eta.high.zero[<<idx>>]-<<level>>.eta.low.zero[<<idx>>]
", .open = "<<", .close = ">>", .trim = FALSE)
}

#' Build the exposure-effect prior for one taxonomic level
#'
#' Generates `for(p in 1:P) { <level>.beta[idx,p] ~ dnorm(...) }`: each
#' level's exposure effect is drawn from a normal distribution centered on a
#' value borrowed from the parent level, computed via `inprod()` over the
#' parent's beta and the taxonomy indicator matrix (Equation 6) - both the
#' count-model and zero-inflation components. If `parent` is `NULL` (only
#' phylum, the top of the hierarchy), there's nothing to borrow from, so the
#' prior falls back to a fixed `dnorm(0, tau)` instead.
#'
#' Shared by every level including species, via [.bahzing_level_block()] and
#' [.bahzing_species_block()] respectively - see `mu_owner` below for why
#' species needs one extra parameter to reuse this.
#'
#' @param level Character. This level's variable-name prefix (e.g. `"genus"`,
#'   `"species"`), used for `<level>.beta`/`<level>.tau`.
#' @param idx Character. This level's JAGS loop index variable (e.g. `"g.r"`,
#'   `"r"`).
#' @param parent Character or `NULL`. The parent level's variable-name prefix
#'   to borrow the prior mean from (e.g. `"family"` for genus). `NULL` only
#'   for phylum, which has no parent.
#' @param parent_R Character or `NULL`. The JAGS scalar holding the parent
#'   level's taxon count (e.g. `"Family.R"`), used to size the `inprod()`.
#'   Required whenever `parent` is supplied.
#' @param parent_data Character or `NULL`. The taxonomy indicator matrix
#'   linking this level's taxa to the parent's (e.g. `"FamilyData"`).
#'   Required whenever `parent` is supplied.
#' @param mu_owner Character. Whose name labels the intermediate "borrowed
#'   mean" node (`mu.<mu_owner>`). Defaults to `parent`, matching the
#'   original model's convention for genus/family/order/class (e.g. genus's
#'   node is `mu.family`, named after *its* parent). Species is the one
#'   exception in the original model - its node is `mu.species` (named after
#'   *itself*, not `mu.genus`) - so [.bahzing_species_block()] overrides this
#'   explicitly to `"species"`. This is a naming inconsistency inherited from
#'   the original model, preserved rather than "fixed", since `mu.*` is never
#'   an monitored/extracted node - only its label differs, not the value.
#' @param comment Character or `NULL`. If supplied, prepended as a standalone
#'   comment line above the `for(p in 1:P)` loop (species has one in the
#'   original, `"prior on exposure effects"`; genus/family/order/class/phylum
#'   don't have one at all).
#' @return A length-1 character string ending in a trailing newline.
#' @keywords internal
#' @noRd
.bahzing_exposure_prior_block <- function(level, idx, parent = NULL, parent_R = NULL,
                                            parent_data = NULL, mu_owner = parent,
                                            comment = NULL) {
  # plain paste0, not glue - glue's default .trim strips leading whitespace
  # and the trailing newline, which would collapse this onto the same line
  # as whatever follows it.
  comment_line <- if (!is.null(comment)) paste0("      # ", comment, "\n") else ""
  body <- if (!is.null(parent)) {
    glue::glue(
"        <<level>>.beta[<<idx>>,p] ~ dnorm(mu.<<mu_owner>>[<<idx>>,p],<<level>>.tau[<<idx>>])
        mu.<<mu_owner>>[<<idx>>,p] <- inprod(<<parent>>.beta[1:<<parent_R>>,p], <<parent_data>>[<<idx>>,1:<<parent_R>>])
        #Zero inflation component
        <<level>>.beta.zero[<<idx>>,p] ~ dnorm(mu.<<mu_owner>>.zero[<<idx>>,p],<<level>>.tau.zero[<<idx>>])
        mu.<<mu_owner>>.zero[<<idx>>,p] <- inprod(<<parent>>.beta.zero[1:<<parent_R>>,p], <<parent_data>>[<<idx>>,1:<<parent_R>>])
", .open = "<<", .close = ">>", .trim = FALSE)
  } else {
    glue::glue(
"        <<level>>.beta[<<idx>>,p] ~ dnorm(0, <<level>>.tau[<<idx>>])
        #Zero inflation component
        <<level>>.beta.zero[<<idx>>,p] ~ dnorm(0, <<level>>.tau.zero[<<idx>>])
", .open = "<<", .close = ">>", .trim = FALSE)
  }
  glue::glue(
"<<comment_line>>      for(p in 1:P) {
<<body>>      }
", .open = "<<", .close = ">>", .trim = FALSE)
}

#' Build one taxonomic level's full JAGS block (genus/family/order/class/phylum)
#'
#' Assembles one level's complete `for(idx in 1:level_R) { ... }` block:
#' a comment header, the exposure-effect prior
#' ([.bahzing_exposure_prior_block()]), and the precision-prior/g-estimation
#' lines ([.bahzing_precision_gestimation_block()]). Handles both
#' genus/family/order/class (`parent` supplied - borrows its prior mean from
#' the parent level) and phylum (`parent = NULL` - the top of the hierarchy,
#' with nothing to borrow from). These two cases used to be separate
#' near-duplicate functions (`.bahzing_mid_level_block()` /
#' `.bahzing_top_level_block()`) before being unified here.
#'
#' Not used for species - species has its own likelihood/dispersion/
#' intercept/covariate structure with no equivalent at any other level, so it
#' gets its own function, [.bahzing_species_block()], which reuses the same
#' two shared sub-block generators this function uses.
#'
#' @param level Character. This level's variable-name prefix, e.g. `"genus"`.
#' @param level_label Character. Human-readable label for this level's
#'   comment header, e.g. `"Genus"` (produces `# Genus level`).
#' @param idx Character. This level's JAGS loop index variable, e.g. `"g.r"`.
#' @param level_R Character. The JAGS scalar holding this level's taxon
#'   count, e.g. `"Genus.R"` - bounds the `for(idx in 1:level_R)` loop.
#' @param parent Character or `NULL`. The parent level's variable-name prefix
#'   (e.g. `"family"` for genus). `NULL` only for phylum.
#' @param parent_R Character or `NULL`. The parent level's taxon-count
#'   scalar (e.g. `"Family.R"`). Required whenever `parent` is supplied.
#' @param parent_data Character or `NULL`. The taxonomy indicator matrix
#'   linking this level to the parent (e.g. `"FamilyData"`). Required
#'   whenever `parent` is supplied.
#' @return A length-1 character string: this level's complete JAGS block.
#' @keywords internal
#' @noRd
.bahzing_level_block <- function(level, level_label, idx, level_R,
                                  parent = NULL, parent_R = NULL, parent_data = NULL) {
  exposure_prior  <- .bahzing_exposure_prior_block(level, idx, parent, parent_R, parent_data)
  precision_gestim <- .bahzing_precision_gestimation_block(level, idx)

  glue::glue(
"    # <<level_label>> level
    for(<<idx>> in 1:<<level_R>>) {
<<exposure_prior>><<precision_gestim>>    }

", .open = "<<", .close = ">>", .trim = FALSE)
}

#' Build the species-level JAGS block - the base/likelihood level
#'
#' Species is the one level tied directly to the observed data, so unlike
#' [.bahzing_level_block()] (genus/family/order/class/phylum, which are pure
#' prior/hierarchy blocks), this one also includes: the zero-inflated
#' negative binomial likelihood (`Y ~ dnegbin(...)`, `zero ~ dbern(...)`),
#' the dispersion prior (`disp`), the intercept priors (`alpha`/
#' `alpha.zero`), and - when `has_covar` is `TRUE` - the covariate terms and
#' their priors (`delta`/`delta.zero`). This is what collapses the original
#' model's separate with/without-covariates literal text into one function.
#'
#' Reuses the same two shared sub-block generators every other level uses:
#' [.bahzing_exposure_prior_block()] (with `mu_owner = "species"` and
#' `comment = "prior on exposure effects"`, matching the original model's
#' species-specific naming/comment conventions) and
#' [.bahzing_precision_gestimation_block()].
#'
#' @param has_covar Logical. Whether covariates are included in the model.
#'   When `TRUE`, adds the `delta`/`delta.zero` covariate terms to the
#'   likelihood and their `dnorm(0, 1.0E-02)` priors; when `FALSE`, omits
#'   them entirely (matching the original model's with/without-covariates
#'   variants).
#' @return A length-1 character string: the complete species-level JAGS
#'   block.
#' @keywords internal
#' @noRd
.bahzing_species_block <- function(has_covar) {
  covar_lambda <- if (has_covar) " + inprod(delta[r, 1:Q], W[i,1:Q])" else ""
  covar_pi     <- if (has_covar) " + inprod(delta.zero[r, 1:Q], W[i,1:Q])" else ""
  covar_prior  <- if (has_covar)
"
      # prior on covariate effects
      for(q in 1:Q) {
        delta[r,q] ~ dnorm(0, 1.0E-02)
        delta.zero[r,q] ~ dnorm(0, 1.0E-02)
      }
" else ""

  exposure_prior <- .bahzing_exposure_prior_block("species", "r", parent = "genus",
                                                    parent_R = "Genus.R", parent_data = "GenusData",
                                                    mu_owner = "species",
                                                    comment = "prior on exposure effects")
  precision_gestim <- .bahzing_precision_gestimation_block("species", "r")

  glue::glue(
"    for(r in 1:R) {
      for(i in 1:N) {
        Y[i,r] ~ dnegbin(mu[i,r], disp[r])
        mu[i,r] <- disp[r]/(disp[r]+(1-zero[i,r])*lambda[i,r]) - 0.000001*zero[i,r]
        log(lambda[i,r]) <- alpha[r] + inprod(species.beta[r,1:P], X.q[i,1:P])<<covar_lambda>>

        # zero-inflation
        zero[i,r] ~ dbern(pi[i,r])
        logit(pi[i,r]) <- alpha.zero[r] + inprod(species.beta.zero[r,1:P], X.q[i,1:P])<<covar_pi>>
      }
      # prior on dispersion parameter
      disp[r] ~ dunif(0,50)

      # prior on intercept
      alpha[r] ~ dnorm(0, 1.0E-02)
      alpha.zero[r] ~ dnorm(0, 1.0E-02)
<<covar_prior>>
<<exposure_prior>>
<<precision_gestim>>    }

", .open = "<<", .close = ">>", .trim = FALSE)
}

#' Assemble the complete JAGS model text for BaHZING_Model()
#'
#' The top-level entry point for the templating system: builds all six
#' taxonomic-level blocks (species via [.bahzing_species_block()]; genus,
#' family, order, class, and phylum via [.bahzing_level_block()], each
#' passed its own parent level's info) and concatenates them into one
#' complete `model { ... }` string, ready to pass to
#' `jags.model(file = textConnection(...))`. Called once per
#' `BaHZING_Model()` invocation as
#' `.bahzing_build_model_text(has_covar = !is.null(covar))`.
#'
#' @param has_covar Logical. Whether the model includes covariates - passed
#'   straight through to [.bahzing_species_block()], the only level that
#'   varies its text based on this (see that function for why).
#' @return A length-1 character string: the complete JAGS model text, from
#'   `model {` through the final closing `}`.
#' @keywords internal
#' @noRd
.bahzing_build_model_text <- function(has_covar) {
  species <- .bahzing_species_block(has_covar)
  genus   <- .bahzing_level_block("genus",  "Genus",  "g.r", "Genus.R",  "family", "Family.R", "FamilyData")
  family  <- .bahzing_level_block("family", "Family", "f.r", "Family.R", "order",  "Order.R",  "OrderData")
  order   <- .bahzing_level_block("order",  "Order",  "o.r", "Order.R",  "class",  "Class.R",  "ClassData")
  class   <- .bahzing_level_block("class",  "Class",  "c.r", "Class.R",  "phylum", "Phylum.R", "PhylumData")
  phylum  <- .bahzing_level_block("phylum", "Phylum", "p.r", "Phylum.R")
  paste0("model {\n", species, genus, family, order, class, phylum, "  }")
}

BaHZING_Model <- function(formatted_data,
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
                          ROPE_range = c(-0.1, 0.1),
                          seed = NULL,
                          return_posterior = FALSE) {

  # JAGS has its own RNG, independent of R's - set.seed() alone has no effect
  # on it. Reproducible runs require explicit .RNG.name/.RNG.seed per chain
  # passed as `inits` to jags.model().
  jags_inits <- if (!is.null(seed)) {
    lapply(seq_len(n.chains), function(i) {
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed + i)
    })
  } else {
    NULL
  }

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
  if(exposure_standardization=="quantile" &
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
  # Matches both the greengenes-style "k__" (kingdom) and SILVA/QIIME2-style
  # "d__" (domain) taxonomy prefixes, since input data may use either.
  Y <- exposure_covar_dat[, grep("^[kd]__", names(exposure_covar_dat))]
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
  L <- exposure_covar_dat[, grep("^[kd]__", names(exposure_covar_dat))]
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

  # 6. Run Model ----
  BHRM.microbiome <- .bahzing_build_model_text(has_covar = !is.null(covar))

  jdata <- list(N=N, Y=Y, R=R, X.q=X.q, P=P,
                GenusData=GenusData, Genus.R=Genus.R,
                Family.R=Family.R, FamilyData=FamilyData,
                Order.R=Order.R, OrderData=OrderData,
                Class.R=Class.R, ClassData=ClassData,
                Phylum.R=Phylum.R, PhylumData=PhylumData,
                profiles=profiles, L=L)
  if (!is.null(covar)) {
    jdata$Q <- Q
    jdata$W <- W
  }

  var.s <- c("species.beta", "genus.beta", "family.beta", "order.beta",
             "class.beta", "phylum.beta", "species.beta.zero",
             "genus.beta.zero", "family.beta.zero", "order.beta.zero",
             "class.beta.zero", "phylum.beta.zero","species.psi","genus.psi",
             "family.psi","order.psi","class.psi","phylum.psi",
             "species.psi.zero","genus.psi.zero","family.psi.zero",
             "order.psi.zero","class.psi.zero","phylum.psi.zero",
             "disp")

  ### Run JAGs Estimation, in parallel across the n.chains independent chains ----
  # The n.chains chains are statistically independent - previously run
  # sequentially inside one jags.model(..., n.chains=n.chains) call on a
  # single core. Each chain is now compiled and sampled in its own forked
  # process via parallel::mclapply, one core per chain.
  if (is.null(jags_inits)) {
    # Without an explicit seed, each chain still needs a distinct RNG seed
    # decided here in the parent process before forking - a forked child's
    # RNG state is otherwise copied identically from the parent, which would
    # give every "independent" chain the same draws.
    chain_seeds <- sample.int(1e6, n.chains)
    jags_inits <- lapply(chain_seeds, function(s) {
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = s)
    })
  }

  run_one_chain <- function(i) {
    m <- jags.model(file=textConnection(BHRM.microbiome), data=jdata,
                    inits=list(jags_inits[[i]]), n.chains=1, n.adapt=n.adapt,
                    quiet=TRUE)
    update(m, n.iter=n.iter.burnin, progress.bar="none")
    chain <- coda.samples(model=m, variable.names=var.s, n.iter=n.iter.sample,
                          thin=1, progress.bar="none")
    chain[[1]]
  }

  chain_list <- parallel::mclapply(seq_len(n.chains), run_one_chain,
                                   mc.cores = min(n.chains, parallel::detectCores()))
  # mclapply returns a try-error object per-fork on failure instead of
  # propagating it - surface that clearly instead of failing downstream with
  # a confusing coda/as.mcmc.list error. A try-error is a character vector
  # (not a condition object), with the actual condition stashed in its
  # "condition" attribute, so conditionMessage() can't be called on it
  # directly - that itself errors ("$ operator is invalid for atomic
  # vectors"), masking the real failure.
  failed <- vapply(chain_list, function(x) inherits(x, "try-error"), logical(1))
  if (any(failed)) {
    msgs <- vapply(chain_list[failed], function(x) {
      cond <- attr(x, "condition")
      if (!is.null(cond)) conditionMessage(cond) else as.character(x)
    }, character(1))
    stop("Parallel chain(s) failed: ", paste(unique(msgs), collapse = "; "))
  }
  model.fit <- coda::as.mcmc.list(chain_list)

  # 7. summarize results -------------------------------------------------------
  ## Calculate Mean, SD, and quantiles ----
  r <- summary(model.fit)
  results <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

  ## Calculate HPD intervals (not included in current version) ------
  # x1 <- HPDinterval(model.fit[[1]], prob = 0.95)  %>% as.data.frame()

  ## Calculate p-values  ------
  # Operate on the native matrix (avoids coercing a large posterior matrix to
  # a data frame), and compute p_direction/p_rope/p_map in a single pass per
  # column instead of three separate apply() scans over the same columns.
  post_dist <- as.matrix(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]])), drop = FALSE]

  p_stats <- apply(post_dist, 2, function(x) {
    c(pdir  = as.numeric(p_direction(x, threshold = 0.05)),
      prope = p_rope(x, rope = ROPE_range)$p_ROPE,
      pmap  = as.numeric(p_map(x)))
  })

  # Get dataframe of p-values
  p_value_df <- data.frame(name = colnames(post_dist), t(p_stats), row.names = NULL)

  # # "Significant" example
  # bayestestR::describe_posterior(post_dist$`family.beta[1,1]`)
  # (pdirection <- bayestestR::p_direction(post_dist$`family.beta[1,1]`, as_p = TRUE)) %>% as.numeric()
  # (prope <- bayestestR::p_rope(post_dist$`family.beta[1,1]`, rope = c(-0.1, 0.1), as_p = TRUE) %>% as.numeric())
  # (pmap <- bayestestR::p_map(post_dist$`family.beta[1,1]`)) %>% as.numeric()
  # (ps <- bayestestR::p_significance(post_dist$`family.beta[1,1]`, rope = c(-0.1, 0.1), as_p = TRUE))
  #
  # # "Non-significant" example
  # bayestestR::describe_posterior(post_dist$`species.beta[189,2]`)
  # bayestestR::p_direction(post_dist$`species.beta[189,2]`, as_p = TRUE)
  # bayestestR::p_rope(post_dist$`species.beta[189,2]`, rope = c(-0.1, 0.1))

  ## Create "component" variable ----
  rn_results <- rownames(results)
  results <- results %>%
    mutate(
      component=case_when(
        grepl("zero",rn_results) ~ "Zero-inflation model coefficients",
        grepl("beta",rn_results) ~ "Count model coefficients",
        grepl("psi",rn_results)  ~ "Count model coefficients",
        grepl("disp",rn_results) ~ "Dispersion",
        #grepl("omega",rn_results) ~ "Omega",
        TRUE ~ "Other"))

  # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
  results <- results %>%
    # mutate(bci_inc_zero = case_when(
    #   grepl("disp",rownames(.)) ~ NA_character_,
    #   grepl("omega",rownames(.)) ~ NA_character_,
    #   (X2.5.<0 & X97.5.<0) | (X2.5.>0 & X97.5.>0) ~ "BCI ",
    #   TRUE ~ "N.S.")) %>%
    rename(estimate = Mean,
           bci_lcl = X2.5.,
           bci_ucl = X97.5.)

  # Combine results df with p-value df
  results$name = rownames(results)
  results2 <- left_join(results, p_value_df, by = "name")
  rownames(results2) <- results2$name
  results2 <- results2 %>% select(-name)
  # Rownames are fixed from here through the filter() below (no reordering
  # or row-adding operations occur in between), so cache once instead of
  # recomputing rownames(results2) on every grepl() call below.
  rn_results2 <- rownames(results2)

  # Calculate odds ratios-- removed --
  # results2 <- results2 %>%
  #   mutate(OR=exp(Mean),
  #          OR.ll=exp(X2.5.),
  #          OR.ul=exp(X97.5.))

  #Format output names
  phylum   <- colnames(PhylumData)
  class    <- colnames(ClassData)
  order    <- colnames(OrderData)
  family   <- colnames(FamilyData)
  genus    <- colnames(GenusData)
  species  <- colnames(Y)
  exposure <- colnames(X)
  results2$taxa_index <- str_remove(rn_results2,"..*\\[")
  results2$taxa_index <- str_remove(results2$taxa_index,",.*$")
  results2$taxa_index <- str_remove(results2$taxa_index,"]")
  results2$taxa_index <- as.numeric(results2$taxa_index)
  results2$Exposure.Index <- str_remove(rn_results2,"..*\\,")
  results2 <- results2 %>%
    mutate(Exposure.Index=ifelse(
      grepl("disp",Exposure.Index)|
        grepl("omega",Exposure.Index)|
        grepl("psi",Exposure.Index),
      NA,Exposure.Index))
  results2$Exposure.Index <- str_remove(results2$Exposure.Index,"]")
  results2$Exposure.Index <- as.numeric(results2$Exposure.Index)


  results2 <- results2 %>%
    mutate(taxa_full=case_when(
      grepl("phylum",rn_results2) ~ paste0(phylum[taxa_index]),
      grepl("class"  ,rn_results2) ~ paste0(class[taxa_index]),
      grepl("order"  ,rn_results2) ~ paste0(order[taxa_index]),
      grepl("family" ,rn_results2) ~ paste0(family[taxa_index]),
      grepl("genus"  ,rn_results2) ~ paste0(genus[taxa_index]),
      grepl("species",rn_results2) ~ paste0(species[taxa_index])),
      domain=case_when(
        grepl("phylum" ,rn_results2)  ~ "Phylum",
        grepl("class"  ,rn_results2)   ~ "Class",
        grepl("order"  ,rn_results2)   ~ "Order",
        grepl("family" ,rn_results2)  ~ "Family",
        grepl("genus"  ,rn_results2)   ~ "Genus",
        grepl("species",rn_results2) ~ "Species"),
      exposure=paste0(exposure[Exposure.Index]))

  # Modify Exposure variable
  results2 <- results2 %>%
    mutate(exposure=case_when(
      grepl("psi",rn_results2) ~ "Mixture",
      grepl("disp",rn_results2) ~ "Dispersion",
      #grepl("omega",rn_results2) ~ "Omega",
      TRUE ~ exposure))

  # Get Taxa and domain information for
  results2 <- results2 %>%
    mutate(taxa_full=ifelse(grepl("disp", rn_results2),paste0(species[taxa_index]),taxa_full),
           taxa_full=ifelse(grepl("omega",rn_results2),paste0(species[taxa_index]),taxa_full),
           taxa_name = sub(".*__", "", taxa_full),
           domain = ifelse(exposure == "Dispersion" | exposure == "Omega",
                           "Species", domain))

  # Remove "disp" and "omega" estimates
  if(!return_all_estimates){
    results2 <- results2 %>%
      filter(!grepl("disp",rn_results2),
             !grepl("omega",rn_results2))
  }

  # Remove rownames
  rownames(results2) <- NULL

  # Select final variables
  results2 <- results2 %>%
    select(taxa_full, taxa_name, domain, exposure,component,
           estimate,bci_lcl,bci_ucl,pdir,prope,pmap)

  if (return_posterior) {
    return(list(results = results2, posterior = model.fit))
  }
  return(results2)
}

