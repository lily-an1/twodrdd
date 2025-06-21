

#' Make the sentinels.
#'
#' Generate a grid of sentinels along the treatment boundary (assumed
#' to be where one of ratings is 0) and number the sentinels.
#'
#' @param dat Data to be analyzed.
#' @param n_sentinel The number of sentinels per side.
#' @param fixed_sent TRUE if user wants fixed sentinel locations between 0-4
#'  (approx the width of test sentinel bands in the simulation) on each running variable.
make_sentinels <- function( dat, n_sentinel, fixed_sent = FALSE) {

  rating2 <- NULL

  stopifnot( is.data.frame(dat) )

  sent2 = NA
  sent1 = NA
  if (fixed_sent == TRUE) {
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


#' Sentinel weights
#'
#' Calculate weights of the sentinels.
#' Here we fit a multivariate normal and then do weights as density
#' given the assumed multivariate normal.
#'
#' @param x The location of sentinels, sentinels$rating1.
#' @param y The location of sentinels, sentinels$rating2.
#' @param data The data to be analyzed.
calc_weights <- function(x, y, data, method = c( "mvnorm", "kernel" ) ) {

  method = match.arg( method )

  if ( method == "kernel" ) {
    # 2D kernel density estimate over rating1 and rating2
    kde <- MASS::kde2d(data$rating1, data$rating2, n = 100)

    # Interpolate density values at (x, y)
    dens_vals <- fields::interp.surface(kde, cbind(x, y))

    dens_vals / sum(dens_vals)

  } else {

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
}





#' Drop sentinels with low weight.
#'
#' Given a list of sentinels, drop those with very low weight (using
#' the weight column).  Will drop smallest sentinels in order up to
#' the drop proportion, but never exceed that threshold in total mass
#' dropped.
#'
#' @param GP_res Gaussian process regression result to be passed in.
#' @param drop_prop Proportion of weight under which we drop the sentinel.
#'
#' @export
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
