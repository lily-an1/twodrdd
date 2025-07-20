

#' Gaussian process predictions
#'
#' Innermost function - runs the laGP model to fit gaussian process
#' and uses it to predict for each (x,y) pair of sentinels.
#'
#' @param sentinels Data frame holding the sentinels.
#' @param data The data to be analyzed, passed.
#' @param method What type of GP to call, passed.
#' @param startnum Parameter to Gaussian process estimator.
#' @param endnum Parameter to Gaussian process estimator.
#' @param prefix Used in Sigma column names.
#'
#' @return Data frame with the sentinels and their GP predictions.  If
#'   variance-covariance available, it will add that in as extra
#'   columns with prefix `prefix`.
calc_one_point_GP <- function( sentinels, data,
                               method = c( "new", "aGP", "laGP", "iso" ),
                               startnum = NULL, endnum = NULL,
                               prefix = "var." ) {

  rating1 <- rating2 <- NULL
  method = match.arg( method )
  stopifnot( !is.null( data$Y ) )

  x = sentinels$rating1
  y = sentinels$rating2
  XX = data.frame(x, y)

  X = dplyr::select(data,rating1,rating2)


  GPmodel = NA
  Sigma = NULL
  fit_info = NULL

  if ( method == "new" ) {
    # This fits the actual GP rather than doing an approximation
    # (more computationally expensive).
    da <- laGP::darg(NULL, X)
    ga <- laGP::garg(list(mle=TRUE), data$Y)

    gpi <- laGP::newGPsep(X, data$Y, d=da$start, g=ga$start, dK=TRUE)
    #d = lengthscale, g = nugget
    mle <- laGP::mleGPsep( gpi, param="both",
                           tmin=c(da$min, ga$min),
                           tmax=c(da$max, ga$max),
                           ab=c(da$ab, ga$ab),
                           maxit=1000 )

    GPmodel <- laGP::predGPsep(gpi, XX, lite = FALSE, nonug = TRUE)

    laGP::deleteGPsep(gpi)

    Sigma = GPmodel$Sigma
    GPmodel$se = sqrt( diag( GPmodel$Sigma ) )
    fit_info = mle

  } else if ( method == "iso" ) {

    # Initial values (use your own logic/da/garg if you want)
    d_init <- 1  # starting value for lengthscale
    g_init <- 1e-4 # starting value for nugget

    # Create global "separable" GP with all ds equal
    gpi <- laGP::newGPsep(X, data$Y, d = rep(d_init, ncol(X)), g = g_init, dK = TRUE)

    # -----
    # Custom joint log-likelihood function: params is c(single_d, g)
    loglik_fun <- function(params) {
      dval <- params[1]
      gval <- params[2]
      # Bound checks (optional)
      if (dval <= 0 || gval <= 0) return(-1e20)
      laGP::updateGPsep(gpi, d = rep(dval, ncol(X)), g = gval)
      -laGP::llikGPsep(gpi)  # negative loglik for minimization
    }

    # Optimize jointly:
    res <- optim(par = c(d_init, g_init), fn = loglik_fun,
                 method = "L-BFGS-B", lower = c(1e-6, 1e-8))
    d_hat <- res$par[1]
    g_hat <- res$par[2]

    # Update the model with MLE params
    laGP::updateGPsep(gpi, d = rep(d_hat, ncol(X)), g = g_hat)

    GPmodel <- laGP::predGPsep(gpi, XX, lite = FALSE)
    GPmodel$se <- sqrt(GPmodel$var)

    laGP::deleteGPsep(gpi)

    Sigma = GPmodel$Sigma
    GPmodel$se = sqrt( diag( GPmodel$Sigma ) )
    GPmodel$da1 = d_hat
    GPmodel$da2 = NA


  } else if ( method == "aGP" ) {
    # This is "Localized Approximate GP Regression For Many Predictive Locations"
    # One consequence: This one does not have the Sigma matrix returned.

    GPmodel <- laGP::aGP(X, data$Y, XX, method="nn", verb = 0 )
    #defaults: start = 6, end = 50, d = NULL, g = 1/10000
    GPmodel$se = sqrt(GPmodel$var)

  } else if ( method == "laGP" ) {
    # This is "Localized Approximate GP Prediction At a Single Input Location"
    stopifnot( !is.null( startnum ) )
    stopifnot( !is.null( endnum ) )

    GPmodel <- laGP::laGP(XX, start=startnum, end=endnum, X, Z=data$Y, method="nn",
                          lite=FALSE, d=NULL, g=NULL) #g default, d default
    Sigma = GPmodel$Sigma
    GPmodel$se = sqrt( diag( GPmodel$Sigma ) )

  } else {
    stop( glue::glue( "unrecognized GP method {method}") )
  }

  output <- data.frame(sentNum = sentinels$sentNum,
                       rating1 = x,
                       rating2 = y,
                       Yhat = GPmodel$mean,
                       se = GPmodel$se )

  if ( !is.null( Sigma ) ) {
    colnames(Sigma) = paste0(prefix, sentinels$sentNum)
    output <- dplyr::bind_cols( output, Sigma )
  }

  if ( !is.null( fit_info ) ) {
    attr( output, "fit_info" ) <- fit_info
  }

  return(output)
}





#' Calculate the treatment curve along the given sentinels.
#'
#' This function fits the GP twice, once for each dataset.
#' Middle function - divides data into Treated & Control before Gaussian process prediction.
#'
#' @param sentinels Data frame of the sentinels.
#' @param data Data to be analyzed, passed.
#' @param method What type of GP to call, passed.
#' @param startnum Parameter to Gaussian process estimator.
#' @param endnum Parameter to Gaussian process estimator.
calc_tx_curve <- function(sentinels, data, method, startnum, endnum) {

  T <- Yhat1 <- Yhat0 <- se0 <- se1 <- sentNum <- rating1 <- rating2 <- estimate <- se <- weight <- weight_p <- NULL

  datY0 = dplyr::filter(data, T == 0)
  Y0s = calc_one_point_GP(
    sentinels,
    data = datY0,
    method = method,
    startnum = startnum, endnum = endnum,
    prefix = "var0.")

  datY1 = dplyr::filter(data, T == 1)
  Y1s = calc_one_point_GP(
    sentinels,
    data = datY1,
    method = method,
    startnum = startnum, endnum = endnum,
    prefix = "var1.")

  results = dplyr::inner_join(Y0s,
                              Y1s,
                              by = c("sentNum", "rating1", "rating2"),
                              suffix = c("0", "1"))  %>%
    dplyr::mutate( estimate=Yhat1-Yhat0,
                   se = sqrt(se0^2 + se1^2) )

  results$weight_p = calculate_precision_weights(results)

  results <- results %>%
    dplyr::mutate(weight = sentinels$weight) %>%
    dplyr::relocate( sentNum, rating1, rating2, Yhat0, se0, Yhat1, se1,
                     estimate, se, weight, weight_p )

  fit1 = attr( Y0s, "fit_info" )
  fit2 = attr( Y1s, "fit_info" )
  attr( results, "fit_info" ) <- list( fit1 = fit1, fit2 = fit2 )

  return(results)
}


#' Use Gaussian process regression
#'
#' Given a two-dimensional dataset with two cutpoints, estimate treatment effects
#' along the boundary.
#' This method uses a Gaussian process regression on the treated and
#' control units to estimate a treatment impact along the boundary
#' defined by having a 0 score in at least one of the two ratings.
#' This function will label units as treated or not, based on the
#' (rating1, rating2) values.
#'
#' @param sampdat Data to analyze.
#' @param n_sentinel The number of sentinels per side.
#' @param method What type of GP to call - newGP, aGP, or laGP (Use
#'   'new' generally, although it is more time intensive, others are
#'   approximations).
#' @param startnum Parameter to Gaussian process estimator.
#' @param endnum Parameter to Gaussian process estimator.
#' @param fixed_sent TRUE if user wants to use fixed sentinel locations between 0-4 on each running variable.
#'
#' @export
gaussianp <- function(sampdat, n_sentinel = 20,
                      method = c( "new", "aGP", "laGP", "iso"  ),
                      residualize = FALSE,
                      startnum, endnum, fixed_sent = FALSE) {

  if ( residualize ) {
    M = stats::lm( Y ~ rating1 + rating2, data=sampdat )

    sampdat <- dplyr::mutate( sampdat,
                              Y_old = Y,
                              Y = resid( M ) )
  }

  method = match.arg( method )

  stopifnot( is.data.frame(sampdat) )

  #Take the min. rating value (score) from each row
  r.c<-apply(sampdat[,c("rating1","rating2")], 1, min)

  ww.gp<-999
  # create dataset with outcome, treatment status defined in gen_sim_data.R,
  # what the score was, rating 1, rating 2
  sampdat<-data.frame( Y=sampdat$Y,
                       T=sampdat$T,
                       r.c,
                       rating1=sampdat$rating1,
                       rating2=sampdat$rating2 )

  n = nrow(sampdat)
  sampdat<-sampdat[r.c>-ww.gp & r.c<ww.gp,] #remove very large ratings
  if ( nrow(sampdat) < n ) {
    warning( "Some observations dropped due to extreme rating values." )
  }

  sentinels = make_sentinels(sampdat, n_sentinel = n_sentinel,
                             fixed_sent = fixed_sent )

  #Drop sentinels with too little data to estimate
  sentinels <- sentinels %>%
    dplyr::mutate(weight = calc_weights( sentinels$rating1,
                                         sentinels$rating2,
                                         data = sampdat )) %>%
    drop_low_weight_sentinels()

  results.gp <- calc_tx_curve(sentinels = sentinels,
                              data = sampdat,
                              method = method,
                              startnum=startnum,
                              endnum=endnum)

  results.gp <- results.gp %>%
    dplyr::mutate(n_sentinel = n_sentinel)

  return(results.gp)
}


#' Return ML parameters of the GP model
#'
#' This function extracts the parameters of the Gaussian process,
#' given the result from gaussianp
#'
#' @param gp_res Result from gaussianp.
#' @export
params <- function( gp_res ) {

  inner_params <- function( a ) {
    tibble(
      length1 = a$theta[1],
      length2 = a$theta[2],
      nugget = a$theta[3],
      its = a$its,
      conv = a$conv,
      msg = a$msg
    )
  }

  # Get the parameters of the GP model
  a = attr( gp_res, "fit_info" )
  if ( is.null( a ) ) {
    return( NULL )
  }

  if ( "fit1" %in% names( a ) ) {
    bind_rows( Co=inner_params(a$fit1),
               Tx=inner_params(a$fit2),
               .id="group" )
  } else {
    inner_params( a )
  }
}


#' Retrieve sentinels' covariance matrix
#'
#' Pull the covariance matrix of the sentinels out of the GP result
#' object (treatment or control side)
#'
#' @param GP_res Gaussian process regression result to be passed in.
#' @param tx Needed to identify Sigmas: treated = 1, control = 0.
get_cov_matrix_side <- function( GP_res, tx = NULL ) {

  pfx = "var."
  if ( !is.null( tx ) ) {
    pfx = paste0( "var", tx, "." )
  }

  # Get the covariance matrix of all the sentinels left in GP_res.
  # (We might have columns of covariances with sentinels we have dropped.)
  keep_sent = GP_res$sentNum
  var_name = paste0( pfx, keep_sent )
  Sigma <- GP_res %>%
    dplyr::select( tidyselect::all_of( var_name ) ) %>%
    as.matrix()

  stopifnot( isSymmetric( unname( Sigma ) ) )
  stopifnot( nrow(Sigma) == ncol(Sigma) )

  Sigma
}


#' Add covariance matrices
#'
#' This is Sigma_b|Y (Equation 6 in Rischard et al. (2021).
#'
#' @param GP_res Gaussian process regression result to be passed in.
get_cov_matrix_tx <- function( GP_res ) {
  Sigma_1 = get_cov_matrix_side( GP_res, tx = 1 )
  Sigma_0 = get_cov_matrix_side( GP_res, tx = 0 )

  Sigma_0 + Sigma_1
}




#' Calculate SEs
#'
#' Calculate the standard error for a weighted average of the
#' sentinels given in the Ys object.
#'
#' @param GP_res Gaussian process regression result to be passed in.
#' @param weight Vector of weights.
#' @param tx Treated side = 1, control = 0.
#'
#' @return Standard error of the average outcome of the sentinels,
#'   weighted by `weight`
calc_Ybar_Var <- function( GP_res, tx = 0, weight ) {

  stopifnot( nrow( GP_res ) == length( weight ) )

  # Get the covariance matrix of all the sentinels left in GP_res.
  # (We might have columns of covariances with sentinels we have dropped.)
  keep_sent = GP_res$sentNum
  var_name = paste0( "var", tx, ".", keep_sent )
  Sigma <- GP_res %>%
    dplyr::select( tidyselect::all_of( var_name ) ) %>%
    as.matrix()

  stopifnot( isSymmetric( unname( Sigma ) ) )
  stopifnot( nrow(Sigma) == ncol(Sigma) )

  SE2 = as.numeric( t(weight) %*% Sigma %*% weight )

  stopifnot( length(SE2) == 1 )

  return( SE2 )
}


#' Precision weights
#'
#' This implements equation 13 from Rischard et al. (2021) to get sentinel
#' weights that correspond to the inverse variance weighting.
#'
#' @param GP_res Gaussian process regression result to be passed in.
calculate_precision_weights <- function( GP_res ) {

  Sigma_bY = get_cov_matrix_tx( GP_res )
  Sigma_inv = solve( Sigma_bY )

  K = nrow( GP_res )
  oneK = rep( 1, K )

  wts = as.numeric( oneK %*% Sigma_inv )
  Z = as.numeric( oneK %*% Sigma_inv %*% oneK )
  wts = wts / Z

  wts
}




#' Calculate the BATE SE
#'
#' Use calc_Ybar_Var function on both the treated and control sides.
#' Weight the sentinels on the boundary as given by 'weight'.
#'
#' @param GP_res Gaussian process regression result to be passed in.
#' @param weight Vector of weights.
calc_AFE_SE <- function( GP_res, weight ) {

  SE2_Y0 = calc_Ybar_Var( GP_res, tx = 0, weight )
  SE2_Y1 = calc_Ybar_Var( GP_res, tx = 1, weight )
  return( sqrt( SE2_Y0 + SE2_Y1 ) )

}


#' Calculate BATE
#'
#' Average the sentinels to get an overall estimated impact along the
#' boundary.
#
#' Produce two estimates: one uses the sentinel weights.  The other
#' uses precision weighting by weighting by the inverse of the standard
#' error (if the variance-covariance of the sentinel estimates is
#' provided).
#'
#' @param GP_res Gaussian process regression result to be passed in.
#' @param calc_SE Will be TRUE for Gaussian process regression and FALSE for loess.
#'
#' @export
calculate_average_impact <- function( GP_res, calc_SE = TRUE ) {

  estimate <- weight <- weight_p <- AFE_wt <- AFE_prec <- AFE <- SE <- parameter <- n_sent_used <- sampsize <- SE_wt <- SE_prec <- NULL

  # If precision weighting is not stored, try to generate some on
  # the fly. (This will be for loess, in general.)
  if ( !("weight_p" %in% names(GP_res)) ) {

    # Make sure we have the variance-covariance matrix
    smpname = paste0( "var1.", GP_res$sentNum )
    if ( all(smpname %in% colnames(GP_res) ) ) {
      GP_res$weight_p = calculate_precision_weights( GP_res )
    } else {
      # We don't--We can make adhoc precision
      # weights without taking correlation into account
      GP_res <- GP_res %>%
        dplyr::mutate(weight_p = 1 / se^2,
                      weight_p = weight_p / sum(weight_p, na.rm=TRUE) )
    }
  }

  if ( calc_SE ) {
    results <- GP_res %>%
      dplyr::summarize( AFE_wt = stats::weighted.mean(estimate, w=weight, na.rm=TRUE),
                        SE_wt = calc_AFE_SE( GP_res, weight ),
                        AFE_prec = stats::weighted.mean( estimate, w = weight_p, na.rm=TRUE ),
                        SE_prec = calc_AFE_SE( GP_res, weight_p ),
                        n_sent_used = dplyr::n(),
                        sampsize = NA,
                        n_sentinel = mean(n_sentinel))

  } else {
    # We might want to skip SE calculation, for example with loess
    results <- GP_res %>%
      dplyr::summarize( AFE_wt = stats::weighted.mean(estimate, w=weight, na.rm=TRUE),
                        SE_wt = NA,
                        AFE_prec = stats::weighted.mean( estimate, w = weight_p, na.rm=TRUE ),
                        SE_prec = NA,
                        n_sent_used = dplyr::n(),
                        sampsize = NA,
                        n_sentinel = mean(n_sentinel))

  }

  results_out <- results %>%
    tidyr::pivot_longer(cols = c(AFE_wt, AFE_prec, SE_wt, SE_prec),
                        names_to = c(".value", "parameter"),
                        names_pattern='(.*)_(.*)' ) %>%
    dplyr::rename(estimate=AFE,
                  se=SE) %>%
    dplyr::mutate(parameter= stringr::str_replace(parameter,"wt","AFE_wt"),
                  parameter= stringr::str_replace(parameter,"prec", "AFE_prec"))

  if ( !calc_SE ) {
    results_out$se = NA
  }

  results_out %>%
    dplyr::relocate( n_sent_used, sampsize, .after = tidyselect::last_col() )
}



