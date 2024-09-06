#from Kristin Porter's "gen.2RRDD_2013_0427.R" code shared Feb 2022
#Modified by Lily An



#' function gen.2RRDD
#'
#' Generate fake data following general structure of Porter et al. (2017), as
#' adjusted by An et al. (2024)
#'
#' @param n = sample size
#' @param s = number of samples
#' @param mu1 = mean of rating 1
#' @param mu2 = mean of rating 2
#' @param sigma.r1 = sd of true rating 1
#' @param sigma.r2 = sd of true rating 2
#' @param rho.r1r2 = correlation between true ratings
#' @param cut1.quantile = quantile for cutting rating 1
#' @param cut2.quantile = quantile for cutting rating 2
#' @param cut1.value = value for cutting rating 1 (if no quantile
#'   passed)
#' @param cut2.value = value for cutting rating 2 (if no quantile
#'   passed).
#' @param param.Y0 = parameters for fn relating ratings and Y0
#' @param param.eff.T = parameters for treatment impact as function of
#'   ratings for full model with terms r1,r2,r1*r2,r1^2,r2^2
#' @param sigma.noise The sd of error term in surface response model
#' @param sigma.E The sd of random error in effects (assumed same for
#'   all params)
#' @param	data.type = "full" or "observed" (ID, rating1, rating2, T, Y)
#'
#' @return: simdat = either full or observed dataset and some
#' descriptive info such as parameter values for observed data.
#'
#' @import mvtnorm
#'
gen.2RRDD <- function(n,
                      s = 1,
                      mu1 = 0,
                      mu2 = 0,
                      sigma.r1 = 1,
                      sigma.r2 = sigma.r1,
                      rho.r1r2 = 0,
                      cut1.quantile = NULL,
                      cut2.quantile = NULL,
                      cut1.value = NULL,
                      cut2.value = NULL,
                      param.Y0 = c( 0, 0, 0, 0, 0, 0 ),
                      tx.func = function( r1, r2 ) { 0.4 },
                      sigma.E = 0,
                      sigma.noise = 1,
                      data.type = "observed",
                      seed = NULL ) {

  if ( !is.null( seed ) ) {
    set.seed(seed)
  }

  # sample id
  sid <- rep(1:s, each = n)

  # generate 2 ratings as a bivariate normal distribution
  cov.r1r2 <- rho.r1r2 * sigma.r1 * sigma.r2
  Sigma = matrix(c(sigma.r1 ^ 2, cov.r1r2, cov.r1r2, sigma.r2 ^ 2), c(2, 2))
  r.true <- MASS::mvrnorm(n * s, mu = c(mu1, mu2), Sigma = Sigma)
  r1.obs <- r.true[, 1]
  r2.obs <- r.true[, 2]


  # assign treatment based on specified cut-point in the observed ratings
  stopifnot(is.null(cut1.quantile) + is.null(cut1.value) == 1)
  stopifnot(is.null(cut2.quantile) + is.null(cut2.value) == 1)

  # in sim, using cutoff location on the N(0,1) distribution
  if (is.null(cut1.quantile) == FALSE) {
    stopifnot( cut1.quantile > 0 && cut1.quantile < 1 )
    cut1 <- qnorm( cut1.quantile, mean = mu1, sd = sigma.r1 )
  } else {
    stopifnot( cut2.quantile > 0 && cut2.quantile < 1 )
    stopifnot( !is.null(cut1.value ) )
    cut1 <- cut1.value
  }
  if (is.null(cut2.quantile) == FALSE) {
    stopifnot( cut2.quantile > 0 && cut2.quantile < 1 )
    cut2 <- qnorm( cut2.quantile, mean = mu2, sd = sigma.r2 )
  } else {
    stopifnot( cut1.quantile > 0 && cut1.quantile < 1 )
    stopifnot( !is.null( cut2.value ) )
    cut2 <- cut2.value
  }

  # Porter - define control as top right quadrant
  # r1fail <- 1 * (r1.obs < cut1)
  # r2fail <- 1 * (r2.obs < cut2)
  # T1 <- 1 * (r1fail == 1 & r2fail == 0)	# upper left quad
  # T3 <- 1 * (r2fail == 1 & r1fail == 0)	# lower rt quad
  # T2 <- 1 * (r1fail == 1 & r2fail == 1) # lower left quad
  # T <- 1 * (T1 == 1 | T2 == 1 | T3 == 1) # treated = 1

  # An - define treatment as top right quadrant
  r1pass <- 1 * (r1.obs>=cut1)
  r2pass <- 1 * (r2.obs>=cut2)
  T1 <- 1 * (r1pass == 1)
  T3 <- 1 * (r2pass == 1)
  T <- 1 * (T1 == 1 & T3 == 1)

  # center ratings around cut
  r1 <- r1.obs - cut1
  r2 <- r2.obs - cut2

  cut1.old <- cut1
  cut2.old <- cut2
  cut1 <- 0
  cut2 <- 0

  # generate potential outcomes
  dm <-
    data.frame(
      intercept = 1,
      r1 = r1,
      r2 = r2,
      r1r2 = r1 * r2,
      r1.2 = r1 ^ 2,
      r2.2 = r2 ^ 2
    )
  # observed ratings and counterfactual:
  #
  # This does a matrix multiply of the design matrix dm (a 6 x n
  # matrix) and the parameter vector param.Y0 (1 x 6)
  Y0 <- apply(sweep(dm, MARGIN = 2, param.Y0, '*'), 1, sum) + rnorm(n, 0, sigma.noise)

  # generate impacts
  eff.err <- rnorm(n * s, 0, sigma.E)
  E.true.T <- tx.func(r1, r2) + eff.err

  # generate potential outcomes under treatment
  Y1 <- Y0 + E.true.T

  # assign potential outcome based on treatment status
  Y.obs <- Y0
  Y.obs[T == 1] <- Y1[T == 1]

  # create dataframe (either full and observed)
  simdat = NULL
  if (data.type == "full") {
    simdat <- data.frame(
      sampleid = sid,
      rating1 = r1,
      rating2 = r2,
      r1pass,
      r2pass,
      cut1,
      cut2,
      T,
      Y0,
      Y1,
      Y.obs
    )
  } else {
    simdat <-
      data.frame(
        sampleID = sid,
        rating1 = r1,
        rating2 = r2,
        T = T,
        Y = Y.obs
      )
  }


  # parameter values
  # sigma
  obs.sigma.r1 <- sd(r1.obs)
  obs.sigma.r2 <- sd(r2.obs)
  obs.sigma.Y <- sd(Y.obs)

  # rho
  obs.rho.r1r2 <- cor(r1.obs, r2.obs)
  obs.rho.r1Y0 <- cor(r1.obs, Y0)
  obs.rho.r2Y0 <- cor(r2.obs, Y0)

  # treatment effects

  # parameter values
  effects = calc_true_effects( r1.obs, r2.obs, E.true.T,
                               rMean = c( mu1 - cut1.old, mu2 - cut2.old ),
                               rSigma = Sigma,
                               tx.func = tx.func )

  # collect all parameter values for outputting
  observed <-
    data.frame(
      cut1,
      cut2,
      obs.sigma.r1,
      obs.sigma.r2,
      obs.sigma.Y,
      obs.rho.r1r2,
      obs.rho.r1Y0,
      obs.rho.r2Y0
    )

  rating = list( cut1 = cut1.old, cut2 = cut2.old,
                 mean = c( mu1, mu2 ),
                 Sigma = Sigma,
                 n = n, s = 1 )

  parameters <- list(effects = effects, observed = observed, rating = rating)


  # stuff to be returned
  out <- list( data = simdat,
               parameters = parameters,
               tx.func = tx.func )

  return(out)

} #end function gen.2RRDD


#' Calculate true effects assuming the running variable distribution
#' is bivariate normal with passed mean and covariance.
#' I.e., (r1, r2) ~ MVN( rMean, rSigma )
calc_true_effects <- function( rating1, rating2,
                               E.true.T,
                               rMean, rSigma,
                               tx.func ) {

  tau.T <- mean(E.true.T) # marginal for whole pop
  tau.T.r2pass <-
    mean(E.true.T[rating2 > 0]) # marginal for those passing r2
  tau.T.r1pass <-
    mean(E.true.T[rating1 > 0]) # marginal for those passing r1

  # Calculate impacts at frontiers assuming bivariate normal density...

  # integrate conditional effect * conditional density
  num.integrand.r1 <- function(xx) {
    tx.func(xx, 0) * dmvnorm(cbind(xx, 0),
                             mean = rMean,
                             sigma = rSigma)
  }
  num.r1 <- integrate(num.integrand.r1, lower = 0, upper = Inf)
  num.integrand.r2 <- function(xx) {
    tx.func(0, xx) * dmvnorm(cbind(0, xx),
                             mean = rMean,
                             sigma = rSigma)
  }
  num.r2 <- integrate(num.integrand.r2, lower = 0, upper = Inf)

  # integrate conditional density
  den.integrand.r1 <- function(xx) {
    dmvnorm(cbind(xx, 0),
            mean = rMean,
            sigma = rSigma)
  }
  den.r1 <- integrate(den.integrand.r1, lower = 0, upper = Inf)
  den.integrand.r2 <- function(xx) {
    dmvnorm(cbind(0, xx),
            mean = rMean,
            sigma = rSigma)
  }
  den.r2 <- integrate(den.integrand.r2, lower = 0, upper = Inf)

  # mean effect along r1 boundary
  tau.r1 <- num.r2$value / den.r2$value
  # mean effect along r2 boundary
  tau.r2 <- num.r1$value / den.r1$value

  # weights for effect among those at either r1 or r2
  wt.denom <- den.r1$value + den.r2$value
  wt.r1 <- den.r2$value / wt.denom
  wt.r2 <- den.r1$value / wt.denom
  tau.T.b <- wt.r1 * tau.r1 + wt.r2 * tau.r2

  effects <- c( ATE = tau.T,
                ATE_r1 = tau.T.r1pass, ATE_r2 = tau.T.r2pass,
                AFE_r1 = tau.r1, AFE_r2 = tau.r2, AFE = tau.T.b,
                wt.r1 = wt.r1, wt.r2 = wt.r2 )

}


#' Given the true mean and covariance, calculate weights of a passed
#' sequence of (r1,r2) pairs.
#'
#' @return List of weights, one for each sentinel.  Sums to 1.
calc_true_sentinel_weights <- function( rating1, rating2,
                                        rMean, rSigma ) {

  wt = dmvnorm( cbind( rating1, rating2 ),
                mean = rMean, sigma = rSigma )

  wt = length(wt) * wt / sum(wt)
}




#' This function generates multiple datasets at once, given the
#' simulation type, cut scores, sample size.
#'
#' Passed parameters allow for changing the data generating process.
#'
#' @param sim ID for the simulation type, following Porter et al.
#' @param n Number of samples per set
#' @param s Number of datasets
#'
#' @export
gen_dat_sim <- function(sim = 1,
                        rho = 0,
                        cut1 = 0.5,
                        cut2 = 0.5,
                        seed = NULL,
                        n = 5000,
                        s = 500) {

  rho.r1r2 <- rho # = correlation between true ratings
  cut1.quantile <- cut1
  # = quantile for cutting rating 1
  cut2.quantile <- cut2
  # = quantile for cutting rating 2

  mu1 <- 0 # = mean of rating 1
  mu2 <- 0 #= mean of rating 2
  sigma.r1 <- 1 #= sd of true rating 1
  sigma.r2 <- 1 #= sd of true rating 2 (always = ratings)
  cut1.value <- NULL # = value for cutting rating 1
  cut2.value <- NULL # = value for cutting rating 2
  sigma.noise <- 1
  sigma.E <- 0 # = sd of random error in individual impacts

  if (sim == 1) {
    # FIRST VALUE ALWAYS 0, COEFF ON R1, COEFF ON R2, ALWAYS 0, COEFF ON R1sq, COEFF ON R2sq
    param.Y0 <- c(0, 0.5, 1, 0, 0, 0)
    param.eff.T <- c(0.4, 0, 0, 0, 0, 0)
  }

  if (sim == 2) {
    param.Y0 <- c(0, 0.5, 1, 0, 2, 1)
    param.eff.T <- c(0.4, 0, 0, 0, 0, 0)
  }

  if (sim == 3) {
    param.Y0 <- c(0, 0.5, 1, 0, 2, 1)
    # COEFF ON T, COEFF OF INTERACTION ON T AND R1, COEFF OF INTERACTION ON T AND R2, ALWAYS 0,
    # COEFF OF INTERACTION ON T AND R1sq, COEFF OF INTERACTION ON T AND R2sq,
    param.eff.T <- c(0.4, -0.1, -0.2, 0, 0, 0)
  }

  if (sim == 4) {
    param.Y0 <- c(0, 0.5, 1, 0, 2, 1)
    # COEFF ON T, COEFF OF INTERACTION ON T AND R1, COEFF OF INTERACTION ON T AND R2, ALWAYS 0,
    # COEFF OF INTERACTION ON T AND R1sq, COEFF OF INTERACTION ON T AND R2sq,
    param.eff.T <- c(0.4, -0.1, -0.2, 0, -0.8, -0.8)
  }


  tx.func = function(r1, r2) {
    pp = param.eff.T
    pp[[1]] + pp[[2]] * r1 + pp[[3]] * r2 + pp[[4]] * r1 * r2 + pp[[5]] *
      r1 ^ 2 + pp[[6]] * r2 ^ 2
  }

  # A new simulation where we have an increasing treatment
  # effect along only one running variable.
  if ( sim == 5 ) {
    param.Y0 <- c(0, 0.5, 1, 0, 2, 1)
    tx.func = function( r1, r2 ) {
      0.2 + r1
    }
  }

  # A simulation that looks like sim 2 but with no treatment effect
  if (sim ==6 ) {
    param.Y0 <- c(0, 0.5, 1, 0, 2, 1)
    tx.func = function( r1, r2 ) {
      0
    }

  }

  gendat <- gen.2RRDD(
    n = n,
    s = s,
    mu1 = mu1,
    mu2 = mu2,
    sigma.r1 = sigma.r1,
    sigma.r2 = sigma.r2,
    sigma.E = sigma.E,
    sigma.noise = sigma.noise,
    rho.r1r2 = rho.r1r2,
    cut1.quantile = cut1.quantile,
    cut2.quantile = cut2.quantile,
    cut1.value = cut1.value,
    cut2.value = cut2.value,
    param.Y0 = param.Y0,
    tx.func = tx.func,
    data.type = "observed",
    seed = seed
  )
}

