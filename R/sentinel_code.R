

#' Make the sentinels.
#'
#' Generate a grid of sentinels along the treatment boundary (assumed
#' to be where one of ratings is 0) and number the sentinels.
#'
#' @param dat Data to be analyzed.
#' @param n_sentinel The number of sentinels per side.
#' @param add_weights Add density weights to the sentinels.
#' @param fixed_sent TRUE if user wants fixed sentinel locations
#'   between 0-max_sent (approx the width of test sentinel bands in the
#'   simulation) on each running variable.
#' @export
make_sentinels <- function( dat, n_sentinel, fixed_sent = FALSE,
                            add_weights = FALSE,
                            drop_low_weight = FALSE,
                            drop_prop = 0.01,
                            max_sent = c( 4, 4 ),
                            ... ) {

  rating2 <- NULL

  stopifnot( is.data.frame(dat) )

  sent2 = NA
  sent1 = NA
  if (fixed_sent == TRUE) {
    sent2 = tibble::tibble(rating1 = 0,
                           rating2 = seq(0, max_sent[[1]],
                                         length.out = n_sentinel)) %>%
      dplyr::arrange(-rating2)
    sent1 = tibble::tibble(rating1 = seq(0, max_sent[[2]],
                                         length.out = n_sentinel),
                           rating2 = 0)

  } else {
    sent2 = tibble::tibble(rating1 = 0,
                           rating2 = seq(0, max(dat$rating2),
                                         length.out = n_sentinel)) %>%
      dplyr::arrange(-rating2)
    sent1 = tibble::tibble(rating1 = seq(0, max(dat$rating1),
                                         length.out = n_sentinel),
                           rating2 = 0)

  }


  sentinels = dplyr::bind_rows(sent2, sent1[-c(1), ])
  sentinels$sentNum = 1:nrow(sentinels)

  if ( add_weights ) {
    sentinels$weight = calc_weights( sentinels$rating1, sentinels$rating2,
                                     data = dat, ... )

    if ( drop_low_weight ) {
      sentinels <- drop_low_weight_sentinels( sentinels, drop_prop = drop_prop )
    }
  }
  sentinels
}


#' Sentinel weights
#'
#' Calculate weights of the sentinels.
#' Here we fit a multivariate normal and then do weights as density
#' given the assumed multivariate normal.
#'
#' @param rating1 The location of sentinels, sentinels$rating1.
#' @param rating2 The location of sentinels, sentinels$rating2.
#' @param data The data to be analyzed.
#'
#' @export
calc_weights <- function(rating1, rating2, data,
                         method = c( "mvnorm", "kernel" ) ) {

  stopifnot( all( (rating1 == 0) | (rating2 == 0) ) )

  method = match.arg( method )

  if ( method == "kernel" ) {
    # 2D kernel density estimate over rating1 and rating2
    kde <- MASS::kde2d(data$rating1, data$rating2, n = 100)

    # Interpolate density values at (rating1, rating2)
    wts <- fields::interp.surface(kde, cbind(rating1, rating2))

  } else {

    mu1 = mean( data$rating1 )
    mu2 = mean( data$rating2 )
    v1 = stats::var( data$rating1 )
    v2 = stats::var( data$rating2 )
    cov12 = stats::cov( data$rating1, data$rating2 )

    wts = mvtnorm::dmvnorm( cbind( rating1, rating2 ),
                            mean = c(mu1, mu2),
                            sigma = matrix( c( v1, cov12, cov12, v2 ), nrow=2 ) )

  }


  # Scale weights to take into account the spacing of the sentinels
  # along each edge. We are assuming the sentinels are equally spaced
  # along the edges, but possibly spaced differently for the different
  # edges.
  dR1 = rating2[[1]] - rating2[[2]]
  stopifnot( all.equal( dR1, rating2[[2]] - rating2[[3]], tolerance = 1e-6 ) )
  n = length(rating1)
  dR2 = rating1[[n]] - rating1[[n-1]]
  stopifnot( all.equal( dR2, rating1[[n-1]] - rating1[[n-2]], tolerance = 1e-6 ) )

  wts[ rating2 == 0 ] = wts[ rating2 == 0 ] * dR2
  wts[ rating1 == 0 ] = wts[ rating1 == 0 ] * dR1
  wts[ rating1 == 0 & rating2 == 0 ] = wts[ rating1 == 0 & rating2 == 0 ] * (dR1+dR2)/( 2 * dR1*dR2 )

  wts <- wts / sum(wts)
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


#' Calculate edge weights.
#'
#' Given a list of sentinels, calculate the total weight on each edge
#' and the total weight overall.
#'
#' @param sents The sentinels to be analyzed.
#' @return A tibble with the total weight, and the weight on each edge
#'
#' @export
calc_edge_weight <- function( sents ) {
  tot <- sum(sents$weight)
  R2 = sum( sents$weight[ sents$rating1 == 0 ] )
  R1 = sum( sents$weight[ sents$rating2 == 0 ] )

  tibble( total = tot, R1 = R1, R2 = R2 )
}


