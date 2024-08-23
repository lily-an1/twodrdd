#
# Given a 2D dataset with two cutpoints, estimate treatment effects
# along the boundary.


#Innermost function - runs the laGP model and prediction on sentinels

#' Generate a grid of sentinels along the treatment boundary (assumed
#' to be where one of ratings is 0)
#'
#' @param dat data to be analyzed
#' @param n_sentinel number of sentinels per side
#' @param fixed_sent fixed sentinel locations between 0-4 on each RV
make_sentinels <- function( dat, n_sentinel, fixed_sent = FALSE) {

    rating2 <- NULL

    stopifnot( is.data.frame(dat) )

    sent2 = NA
    sent1 = NA
    if (fixed_sent == TRUE) {

        #12-07-2022 we decided on 0-4 because that's approx the width of test sentinel bands
        sent2 = tibble::tibble(rating1 = 0, rating2 = seq(0, 4, length.out = n_sentinel)) %>%
            dplyr::arrange(-rating2)
        sent1 = tibble::tibble(rating1 = seq(0, 4, length.out = n_sentinel), rating2 = 0)

    } else {

        sent2 = tibble::tibble(rating1 = 0, rating2 = seq(0, max(dat$rating2), length.out = n_sentinel)) %>%
            dplyr::arrange(-rating2)
        sent1 = tibble::tibble(rating1 = seq(0, max(dat$rating1), length.out = n_sentinel), rating2 = 0)

    }
    sentinels = dplyr::bind_rows(sent2, sent1[-c(1), ])
    sentinels$sentNum = 1:nrow(sentinels)

    sentinels
}


#' Calculate weights of the sentinels.
#'
#' Here we fit a multivariate normal and then do weights as density
#' given the assumed multivariate normal.
#'
#' @param x location of sentinels, sentinels$rating1
#' @param y location of sentinels, sentinels$rating2
#' @param data data to be analyzed
calc_weights <- function(x, y, data) {

    mu1 = mean( data$rating1, )
    mu2 = mean( data$rating2 )
    v1 = stats::var( data$rating1 )
    v2 = stats::var( data$rating2 )
    cov12 = stats::cov( data$rating1, data$rating2 )

    wts = mvtnorm::dmvnorm( cbind( x, y ),
                   mean = c(mu1, mu2),
                   sigma = matrix( c( v1, cov12, cov12, v2 ), nrow=2 ) )

    wts / sum(wts)
}




#' Fit gaussian process and use it to predict for each (x,y) pair
#' passed.
#'
#' @param sentinels Dataframe holding rating1,rating2, which are the
#' @param data data to be analyzed, passed
#' @param method What type of GP to call, passed
#' @param startnum Parameter to gaussian process estimator
#' @param endnum Parameter to gaussian process estimator
#' @param prefix used in Sigma column names
#'   x,y locations of the sentinels.
calc_one_point_GP <- function( sentinels, data,
                               method = c( "new", "aGP", "laGP" ),
                               startnum = NULL, endnum = NULL,
                               prefix = "var." ) {

    rating1 <- rating2 <- NULL
    method = match.arg( method )
    stopifnot( !is.null( data$Y ) )

    x = sentinels$rating1
    y = sentinels$rating2

    X = dplyr::select(data,rating1,rating2) #data is dat.gp by treatment status

    XX = data.frame(x, y) #sentinels
    #safe_aGP <- safely(aGP, otherwise=NA)

    GPmodel = NA
    Sigma = NULL
    if ( method == "new" ) {
        # This actually fits a GP rather than doing an approximation

        da <- laGP::darg(NULL, X)
        ga <- laGP::garg(list(mle=TRUE), data$Y)

        gpi <- laGP::newGPsep(X, data$Y, d=da$start, g=ga$start, dK=TRUE)
        #d = lengthscale, g = nugget
        mle <- laGP::mleGPsep( gpi, param="both", tmin=c(da$min, ga$min),
                         tmax=c(da$max, ga$max), ab=c(da$ab, ga$ab), maxit=1000)

        GPmodel <- laGP::predGPsep(gpi, XX, lite = FALSE, nonug = TRUE)

        laGP::deleteGPsep(gpi)

        Sigma = GPmodel$Sigma
        GPmodel$se = sqrt( diag( GPmodel$Sigma ) )

    } else if ( method == "aGP" ) {
        # This is "Localized Approximate GP Regression For Many
        # Predictive Locations"
        # -- meaning what.
        #
        # One consequence: This one does not have the Sigma matrix
        # returned.

        GPmodel <- laGP::aGP(X, data$Y, XX, method="nn", verb = 0 )
        #defaults: start = 6, end = 50, d = NULL, g = 1/10000
        GPmodel$se = sqrt(GPmodel$var)
        #GPmodel$se = GPmodel$sd/sqrt(length(GPmodel$var))
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

    return(output)
}





#' Calculate the treatment curve along the given sentinels.
#'
#' This function fits the GP twice, once for each dataset.
#'
#' Middle function - divides data into Treated & Control before laGP prediction
#'
#' @param sentinels df of sentinels
#' @param data data to be analyzed, passed
#' @param method What type of GP to call, passed
#' @param startnum Parameter to gaussian process estimator
#' @param endnum Parameter to gaussian process estimator
calc_tx_curve <- function(sentinels, data, method, startnum, endnum) {

    T.c <- Yhat1 <- Yhat0 <- se0 <- se1 <- sentNum <- rating1 <- rating2 <- estimate <- se <- weight <- weight_p <- NULL

    datY0 = dplyr::filter(data, T.c == 0)
    Y0s = calc_one_point_GP(
        sentinels,
        data = datY0,
        method = method,
        startnum = startnum, endnum = endnum,
        prefix = "var0.")

    datY1 = dplyr::filter(data, T.c == 1)
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

    return(results)
}


#' Estimate impacts along boundary given dataset
#'
#' This method uses a gaussian process regression on the treated and
#' control units to estimate a treatment impact along the boundary
#' defined by having a 0 score in at least one of the two ratings.
#'
#' This function will label units as treated or not, based on the
#' (rating1, rating2) values.
#'
#' @param sampdat data to analyze
#' @param n_sentinel number of sentinels per side
#' @param method What type of GP to call - newGP, aGP, or laGP (newGP generally, others are approximations)
#' @param startnum Parameter to gaussian process estimator
#' @param endnum Parameter to gaussian process estimator
#' @param fixed_sent fixed sentinel locations between 0-4 on each RV
gaussianp <- function(sampdat, n_sentinel = 20,
                      method = c( "new", "aGP", "laGP" ),
                      startnum, endnum, fixed_sent = FALSE) {

    method = match.arg( method )

    stopifnot( is.data.frame(sampdat) )

    #Take the min. rating value (score) from each row
    r.c<-apply(sampdat[,c("rating1","rating2")],1,min)

    # If the min. is below 0, mark as treated
    # note: this assumes RVs are centered around their respective cut points.
    T.c<-1*(r.c<0)

    ww.gp<-999
    # create dataset with outcome, "treated" (where treated if score<0),
    # what the score was, rating 1, rating 2
    sampdat<-data.frame( Y=sampdat$Y,
                         T.c,
                         r.c,
                         rating1=sampdat$rating1,
                         rating2=sampdat$rating2 )

    n = nrow(sampdat)
    sampdat<-sampdat[r.c>-ww.gp & r.c<ww.gp,] #remove very large ratings
    if ( nrow(sampdat) < n ) {
        warning( "Some observations dropped due to extreme rating values." )
    }

    sentinels = make_sentinels(sampdat, n_sentinel = n_sentinel, fixed_sent = fixed_sent )

    #1-11-2024 update: Drop sentinels with too little data to estimate
    sentinels <- sentinels %>%
        dplyr::mutate(weight = calc_weights( sentinels$rating1,
                                      sentinels$rating2,
                                      data = sampdat )) %>%
        drop_low_weight_sentinels()

    results.gp <- calc_tx_curve(sentinels, data = sampdat,
                                method = method,
                                startnum=startnum,
                                endnum=endnum)

    return(results.gp)
} # end function gaussianp



#' Given a list of sentinels, drop those with very low weight (using
#' the weight column).  Will drop smallest sentinels up to the drop
#' proportion, but never exceed that threshold.
#'
#' @param GP_res Gaussian process regression result to be passed in
#' @param drop_prop proportion of weight that
drop_low_weight_sentinels <- function( GP_res, drop_prop = 0.01 ) {
    weight <- NULL

    tot = sum( GP_res$weight )
    GP_res$weight = GP_res$weight / sum( GP_res$weight )

    wts = sort( GP_res$weight )
    a = cumsum( wts )
    n_drop = sum( a <= drop_prop )
    n_drop
    wt_cut <- wts[n_drop]

    dplyr::filter( GP_res, weight > wt_cut )
}


#' Pull the covariance matrix of the sentinals out of the GP result
#' object (treatment or control side)
#'
#' @param GP_res Gaussian process regression result to be passed in
#' @param tx treated = 1, control = 0, needed to identify Sigmas
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


#' See maxime's paper for this formula (Equation 6)
#' This is Sigma_b|Y
#'
#' @param GP_res Gaussian process regression result to be passed in
get_cov_matrix_tx <- function( GP_res ) {
    Sigma_1 = get_cov_matrix_side( GP_res, tx = 1 )
    Sigma_0 = get_cov_matrix_side( GP_res, tx = 0 )

    Sigma_0 + Sigma_1
}





#' Calculate the standard error for a weighted average of the
#' sentinels given in the Ys object.
#'
#' @param GP_res Gaussian process regression result to be passed in
#' @param weight vector of weights.
#' @param tx treated side = 1 control = 0
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


#' This implements equation 13 from Maxime paper to get sentinal
#' weights that correspond to the inverse variance weighting.
#'
#' @param GP_res Gaussian process regression result to be passed in
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



#' This implements equation 13 from Maxime paper, and gives the
#' inverse variance weighted BATE with associated estimated standard
#' error.
#'
#' NOTE: The code below does not use this function since it calculates
#' weights via calculate_precision_weights() and then just uses those
#' weights as any other weights for the following calculations.
#'
#' @param GP_res Gaussian process regression result to be passed in
calculate_precision_BATE <- function( GP_res ) {

    Sigma_bY = get_cov_matrix_tx( GP_res )
    Sigma_inv = solve( Sigma_bY )

    K = nrow( GP_res )
    oneK = rep( 1, K )

    wts = as.numeric( oneK %*% Sigma_inv )
    Z = as.numeric( oneK %*% Sigma_inv %*% oneK )
    wts = wts / Z

    BATE_hat = wts %*% t(GP_res$estimate)

    # From Maxime paper, this is a consequence of the usual w'Sigma w
    # formula applied to the precision weighting weights.
    SE2 = 1 / Z

    tibble::tibble( AFE = BATE_hat,
            SE = sqrt( SE2 ),
            n_sent = K )

}



#' Calculate the Standard Error for the Average Frontier Effect
#'
#' Weight the sentinels on the frontier as given by 'weight'.
#'
#' @param GP_res Gaussian process regression result to be passed in
#' @param weight vector of weights
calc_AFE_SE <- function( GP_res, weight ) {
    #Have an SE on the Treatment side, an SE on the Control side
    #sqrt(SE_treat^2 + SE_control^2)

    SE2_Y0 = calc_Ybar_Var( GP_res, tx = 0, weight )
    SE2_Y1 = calc_Ybar_Var( GP_res, tx = 1, weight )
    return( sqrt( SE2_Y0 + SE2_Y1 ) )
}


#' Average the sentinels to get an overall estimated impact along the
#' boundary.
#
#' Produce two estimates: one uses the sentinel weights.  The other
#' uses precision weighting by weighting by the inverse of the standard
#' error (if the variance-covariance of the sentinal estimates is
#' provided).
#'
#' @param GP_res Gaussian process regression result to be passed in
#' @param calc_SE Will be TRUE for GPR and FALSE for loess
calculate_average_impact <- function( GP_res, calc_SE = TRUE ) {

    estimate <- weight <- weight_p <- AFE_wt <- AFE_prec <- AFE <- SE <- parameter <- n_sent <- sampsize <- SE_wt <- SE_prec <- NULL
    #  if (exists("GP_res$error")) {

    # Drop sentinels with too little data to estimate
    # GP_res <- GP_res %>%
    #   dplyr::filter(!is.na(sentNum)) %>%
    #   drop_low_weight_sentinels()

    # If precision weighting is not stored, try to generate some on
    # the fly. (This will be for loess, in general.)
    if ( is.null( GP_res$weight_p ) ) {

        # Make sure we have the variance-covariance matrix
        smpname = paste0( "var1.", GP_res$sentNum )
        if ( all(smpname %in% colnames(GP_res) ) ) {
            GP_res$weight_p = calculate_precision_weights( GP_res )
        } else {
            # We don't--must be loess or something.  We can make adhoc
            # precision weights without taking correlation into
            # account
            GP_res$weight_p = 1 / GP_res$se^2
            GP_res$weight_p = GP_res$weight_p / sum(GP_res$weight_p, na.rm=TRUE)
        }
    }

    if ( calc_SE ) {
        resultskp <- GP_res %>%
            dplyr::summarize( AFE_wt = stats::weighted.mean(estimate, w=weight, na.rm=TRUE),
                              SE_wt = calc_AFE_SE( GP_res, weight ),
                              AFE_prec = stats::weighted.mean( estimate, w = weight_p, na.rm=TRUE ),
                              SE_prec = calc_AFE_SE( GP_res, weight_p ),
                              n_sent = dplyr::n(),
                              sampsize = NA)  #not going to get a sample size for radii along laGP

    } else {
        # We might want to skip SE calculation
        resultskp <- GP_res %>%
            dplyr::summarize( AFE_wt = stats::weighted.mean(estimate, w=weight, na.rm=TRUE),
                              SE_wt = NA,
                              AFE_prec = stats::weighted.mean( estimate, w = weight_p, na.rm=TRUE ),
                              SE_prec = NA,
                              n_sent = dplyr::n(),
                              sampsize = NA)  #not going to get a sample size for radii along laGP

    }

    results_out <- resultskp %>%
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
        dplyr::relocate( n_sent, sampsize, .after = tidyselect::last_col() )
}



