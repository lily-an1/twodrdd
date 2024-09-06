
# Function to estimate treatment curve along a boundary with loess lines
#



#' Given a point (x,y) take all data within radius of the point, and then weight
#' the data in the radius using the triweight function of (1-d^3)^3 (rescaling
#' d, the distance, to have a max of 1 first).  Note the weight is 1 when d=0,
#' and 0 when d=1.
#'
#' Once weights are computed, fit a quad model with a linear interaction (other
#' choices are possible here) to the weighted data.  Then predict for our
#' candidate (x,y) point!
#'
#' Don't fit a model if we have fewer than 6 data points (need sufficient
#' degrees of freedom).
#'
#' @param x one side of sentinels
#' @param y other side of sentinels
#' @param radius radius size around each point
#' @param data data to be analyzed
#' @param min_sample min # of data points within each radius
loess_calc_one_point <- function(x, y, radius, data, min_sample = 8) {

  rating1 <- rating2 <- dist2 <- NULL

  dd = dplyr::mutate(data,
                dist2 = (x - rating1) ^ 2 + (y - rating2) ^ 2) %>%
        dplyr::filter(dist2 <= radius ^ 2)

    if (nrow(dd) < min_sample) {
        return(data.frame(
            Yhat = NA,
            n = nrow(dd),
            eff = 0, # No estimate, no weight.
            se = NA
        ))
    }

    dd <- dd %>%
        dplyr::mutate(dist2 = dist2 / max(dist2),
               w2 = (1 - dist2 ^ (1.5)) ^ 3)

    mod = stats::lm(Y ~ rating1 * rating2, data = dd, weights = dd$w2)
    pd = stats::predict(mod, newdata = data.frame(rating1 = x, rating2 = y), se.fit=T) #, se.fit=TRUE
    eff = ((sum(dd$w2)) ^ 2) / sum(dd$w2 ^ 2)
    data.frame(Yhat = pd$fit,
               n = nrow(dd),
               eff = eff,
               se=pd$se
    )
    }


#' loess_calc_point
#' Run loess_calc_one_point across data
#'
#' @param x one side of sentinels
#' @param y other side of sentinels
#' @param radius radius size around each point
#' @param data data to be analyzed
loess_calc_point <- function(x, y, radius, data) {

    rating1 <- rating2 <- NULL

    pts = purrr::map2_df(x, y, loess_calc_one_point, radius = radius, data = data)
    pts$rating1 = x
    pts$rating2 = y
    pts %>% dplyr::relocate(rating1, rating2)
}



#' Given data and a list of sentinels, calculate predicted Y0 and Y1, and then
#' take the difference as the estimated treatment effect at each sentinel
#'
#' @param sampdat data to be analyzed
#' @param radius radius size around each point
#' @param min_sample min # of data points within each radius
#' @param n_sentinel number of sentinels per side
#'
#' @export
loess2DRDD <- function(sampdat, radius, min_sample, n_sentinel = 20) {

    Yhat1 <- Yhat0 <- se1 <- se0 <- rating1 <- rating2 <- NULL

    sentinels <- make_sentinels( sampdat, n_sentinel = n_sentinel )

    Y0s = loess_calc_point(
        sentinels$rating1,
        sentinels$rating2,
        radius = radius,
        data = dplyr::filter(sampdat, T == 0)
    )
    Y1s = loess_calc_point(
        sentinels$rating1,
        sentinels$rating2,
        radius = radius,
        data = dplyr::filter(sampdat, T == 1)
    )
    Y0s
    Y1s

    Ys = dplyr::inner_join(Y0s,
                    Y1s,
                    by = c("rating1", "rating2"),
                    suffix = c("0", "1"))

    Ys <- Ys %>%
        dplyr::mutate(estimate = Yhat1 - Yhat0,
               se = sqrt(se1^2 + se0^2), #NA, #se = (sqrt(1/sum(1/(se)^2)))
               se_hat = NA, #I don't think we get the covariance matrix from predict.stats::lm
               weight = calc_weights( sentinels$rating1, sentinels$rating2, sampdat ),
               pt = paste0(round(rating1, digits = 1), ",", round(rating2, digits = 1)) )

    Ys$sentNum = 1:nrow(Ys)
    return(Ys)

}

