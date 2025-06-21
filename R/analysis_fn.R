#' Two-dimensional RDD methods
#'
#' Code to analyze single dataset of n observations using our list of
#' estimators. The treated area is the top right quadrant.
#'
#' @param cut1 Cutoff value for running variable 1. Centered at 0 in the simulation's data generating processes.
#' @param cut2 Cutoff value for running variable 2. centered at 0 in the simulation's data generating processes.
#' @param startnum Parameter to the gaussian process estimator.
#' @param endnum Parameter to the gaussian process estimator.
#' @param dat Data to be analyzed.
#' @param n_sentinel The number of sentinels per side.
#' @param include_OLS TRUE to include OLS regression.
#' @param include_BINDING TRUE to include binding score.
#' @param include_FRONTIER TRUE to include pooled frontier.
#'  The results parameter is labeled the AFE. AFE1 and AFE2 are the effects specific to
#'  each running variable boundary.
#' @param include_GP TRUE to include Gaussian process regression.
#' @param include_GP_RES TRUE to include residualized Gaussian process regression.
#' @param include_loess TRUE to include loess.
#' @param include_PARA_LINEAR TRUE to include parametric linear surface response.
#' @param include_PARA_QUAD TRUE to include parametric quadratic surface response.
#'
#' @importFrom stats sd weighted.mean
#' @export
analysis <- function( dat, cut1=0, cut2=0,
                      n_sentinel = 20,
                      startnum=50, endnum=100,
                      include_OLS = TRUE,
                      include_BINDING = FALSE,
                      include_FRONTIER = FALSE,
                      include_GP = FALSE,
                      include_GP_RES = FALSE,
                      include_loess = FALSE,
                      include_PARA_LINEAR = FALSE,
                      include_PARA_QUAD = FALSE) {

    Y <- NULL

    # Center running variables around their cutpoints
    dat$rating1 <- dat$rating1 - cut1
    dat$rating2 <- dat$rating2 - cut2

    ## OLS
    out_ols = tibble::tibble()
    if ( include_OLS ) {
        M = stats::lm( Y ~ rating1 + rating2 + T, data=dat )
        Ms = summary(M)
        out_ols = tibble::tibble( parameter = "AFE",
                          estimate = stats::coef(M)["T"],
                          se = stats::coef(Ms)[4,2],
                          n_sent = NA,
                          n = nrow(dat),
                          sampsize = nrow(dat) )
    }


    ## Binding Score y ~ x1 + x2 | c1 + c2 for a sharp MRDD with two covariates
    out_bs = tibble::tibble()
    if ( include_BINDING ) {
        runfun_bs <- rddapp::mrd_est(data=dat, Y~rating1 + rating2 | rating1 + rating2,
                             method="center",
                             cutpoint=c(cut1, cut2),
                             t.design=c("l","l"))
        tau.afe <- -1*runfun_bs[["center"]][["tau_MRD"]][["est"]][4]
        se.afe <- runfun_bs[["center"]][["tau_MRD"]][["se"]][4]
        n.afe <- runfun_bs[["center"]][["tau_MRD"]][["obs"]][4]
        # [4] is the Optimal bandwidth specification
        out_bs <- c(tau.afe, se.afe, n.afe) # 3 columns - tau, se, n

        both <- rbind(out_bs, out_bs)
        parameter <- c("AFE", "AFE")
        out_bs <- data.frame(parameter, both)
        colnames(out_bs) <- c("parameter", "estimate", "se", "sampsize")
        out_bs <- out_bs %>%
            dplyr::filter(parameter=="AFE")
        out_bs$n = nrow( dat )
    }

    ## Frontier
    out_fr = tibble::tibble()
    if ( include_FRONTIER ) {

        dat$c1 = dat$rating1 > cut1
        dat$c2 = dat$rating2 > cut2

        # Either 0 or all
        if ( sum( dat$c1 ) %in% c(0, nrow(dat) ) ) {
            stop( "c1 is not a good threshold" )
        }
        if ( sum( dat$c2 ) %in% c(0, nrow(dat)) ) {
            stop( "c2 is not a good threshold" )
        }

        tb <- table( c1= dat$c1, c2=dat$c2 )
        tb

        if ( any( tb[,2] <= 10 ) ) {
            # do RDD on c1 only
            runfun_frontier1<-rddapp::rd_est(data=dat, Y~rating1,
                                    cutpoint=cut1,
                                    t.design="l")

            tau.afe1 <- -1*runfun_frontier1[["est"]][4]
            se.afe1 <- runfun_frontier1[["se"]][4]
            n.afe1 <- runfun_frontier1[["obs"]][4]
            # [4] is the Optimal bandwidth specification
            out_fr <- c(tau.afe1, se.afe1, n.afe1) # 3 columns - tau, se, n

            both <- rbind(out_fr, out_fr)
            parameter <- c("AFE1", "AFE1")
            out_fr <- data.frame(parameter, both)
            colnames(out_fr) <- c("parameter", "estimate", "se", "sampsize")
            out_fr <- out_fr %>%
                dplyr::filter(parameter=="AFE1")
            out_fr$n = nrow( dat )

         } else if ( any( tb[2,] <= 10 ) ) {
            # do RDD on c2 only
            runfun_frontier2<-rddapp::rd_est(data=dat, Y~rating2,
                                    cutpoint=cut2,
                                    t.design="l")

            tau.afe2 <- -1*runfun_frontier2[["est"]][4]
            se.afe2 <- runfun_frontier2[["se"]][4]
            n.afe2 <- runfun_frontier2[["obs"]][4]
            # [4] is the Optimal bandwidth specification
            out_fr <- c(tau.afe2, se.afe2, n.afe2) # 3 columns - tau, se, n

            both <- rbind(out_fr, out_fr)
            parameter <- c("AFE2", "AFE2")
            out_fr <- data.frame(parameter, both)
            colnames(out_fr) <- c("parameter", "estimate", "se", "sampsize")
            out_fr <- out_fr %>%
                dplyr::filter(parameter=="AFE2")
            out_fr$n = nrow( dat )

        } else {
            # do RDD on c1 and c2
        runfun_frontier<-rddapp::mrd_est(data=dat, Y~rating1 + rating2 | rating1 + rating2,
                                 method="front",
                                 cutpoint=c(cut1, cut2),
                                 t.design=c("l","l"), boot=50)

        ## AFE using optimal bandwidth (bw selected by CV)
        tau.afe <- -1*runfun_frontier[["front"]][["tau_MRD"]][["est"]][2,3]
        se.afe <- runfun_frontier[["front"]][["tau_MRD"]][["se"]][2,3]
        n.afe <- runfun_frontier[["front"]][["tau_MRD"]][["obs"]][["bw"]][1]

        ## Running variable 1 Optimal bandwidth case (bw selected by CV)
        tau.afe1 <- -1*runfun_frontier[["front"]][["tau_MRD"]][["est"]][2,1]
        se.afe1 <- runfun_frontier[["front"]][["tau_MRD"]][["se"]][2,1]
        n.afe1 <- runfun_frontier[["front"]][["tau_MRD"]][["obs"]][["bw"]][1]

        ## Running variable 2 Optimal bandwidth case (bw selected by CV)
        tau.afe2 <- -1*runfun_frontier[["front"]][["tau_MRD"]][["est"]][2,2]
        se.afe2 <- runfun_frontier[["front"]][["tau_MRD"]][["se"]][2,2]
        n.afe2 <- runfun_frontier[["front"]][["tau_MRD"]][["obs"]][["bw"]][1]

        out1<-c(tau.afe,se.afe,n.afe)
        out2<-c(tau.afe1,se.afe1,n.afe1)
        out3<-c(tau.afe2,se.afe2,n.afe2)

        both <- (rbind(out1, out2, out3))
        parameter<-c("AFE", "AFE1", "AFE2")
        out_fr<-data.frame(parameter, both) %>%
            dplyr::distinct()
        colnames(out_fr)<-c("parameter","estimate","se","sampsize")
        out_fr$n = nrow( dat )
        }
    }

    ## Gaussian process regression
    out_gp = tibble::tibble()
    if ( include_GP ) {
      runfun_gaussianp <- gaussianp(sampdat=dat,
                                    n_sentinel = n_sentinel,
                                    startnum = startnum, endnum=endnum,
                                    method = "new" )

      out_gp = calculate_average_impact( runfun_gaussianp )
      out_gp$n = nrow( dat )

    }

    ## Residualized Gaussian process regression
    out_gp_res = tibble::tibble()
    if (include_GP_RES) {
        M = stats::lm( Y ~ rating1 + rating2, data=dat )

        resdat <- cbind(dat, M$residuals) %>%
            dplyr::mutate(Y_old = Y,
                   Y=M$residuals)

        runfun_gaussianp <- gaussianp(sampdat=resdat,
                                      n_sentinel = n_sentinel,
                                      startnum = startnum, endnum=endnum,
                                      method = "new" )

        out_gp_res = calculate_average_impact( runfun_gaussianp )
        out_gp_res$n = nrow( dat )
    }

    ### Loess approach
    out_loess = tibble::tibble()
    if ( include_loess ) {
      #dat is original data, sampdat is the bootstrapped version
      sampdat <-
        dat %>%
        modelr::bootstrap(50)

      sentinels <- make_sentinels( dat, n_sentinel = n_sentinel )

      sentinels <- sentinels %>%
        dplyr::mutate(weight = calc_weights( sentinels$rating1,
                                      sentinels$rating2,
                                      data = dat )) %>%
        drop_low_weight_sentinels()

      suppressWarnings(out_loess <- purrr::map(sampdat$strap, ~ loess2DRDD(.x,
                                                   radius=mean(sd(dat$rating2), sd(dat$rating1))/2,
                                                   min_sample=min_sample,
                                                   n_sentinel=n_sentinel,
                                                   sentinels=sentinels)))
      out_loess <- dplyr::bind_rows(out_loess, .id = "bootstrap_sample_id")

      out_loess <- out_loess %>%
        dplyr::group_by(bootstrap_sample_id) %>%
        calculate_average_impact(calc_SE=FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(parameter) %>%
        dplyr::summarize(n_sentinel = mean(n_sentinel),
                  se = stats::sd(estimate, na.rm=T),
                  estimate = mean(estimate, na.rm=T),
                  n_sent_used = NA,
                  sampsize = NA)

      out_loess$n = nrow(dat)

    }

    ## Parametric surface response: linear
    out_paraL = tibble::tibble()
    if (include_PARA_LINEAR) {

      sentinels <- make_sentinels( dat, n_sentinel = n_sentinel )
      df <- data.frame(rating1=sentinels$rating1, rating2=sentinels$rating2)

      #Build the model matrix (for predictions, ensures proper variable ordering):
      mm <- stats::model.matrix(~ rating1 + rating2 + rating1:rating2, data=df)

      dat_0 <- dplyr::filter(dat, T==0)
      dat_1 <- dplyr::filter(dat, T==1)

      M0 = estimatr::lm_robust(Y ~ rating1 + rating2 + rating1:rating2, data = dat_0)
      M1 = estimatr::lm_robust(Y ~ rating1 + rating2 + rating1:rating2, data = dat_1)

      #Get estimated coefficients:
      b0 <- stats::coef(M0)
      b1 <- stats::coef(M1)

      #Get predictions for each sentinel point:
      Yhat0 <- as.vector(mm %*% b0)  # Control model
      Yhat1 <- as.vector(mm %*% b1)  # Treated model
      tau_hat <- Yhat1 - Yhat0       # RD effect at each sentinel

      #Extract coefficient covariance matrices:
      vcov0 <- stats::vcov(M0)
      vcov1 <- stats::vcov(M1)

      #vector containing the prediction variance at each sentinel point.
      var0 <- apply(mm, 1, function(x) as.numeric(t(x) %*% vcov0 %*% x))
      var1 <- apply(mm, 1, function(x) as.numeric(t(x) %*% vcov1 %*% x))
      # Variance of treatment effect estimate at each sentinel
      var_tau <- var0 + var1

      #The full covariance matrix of the predictions at all sentinels:
      Sig0 <- mm %*% vcov0 %*% t(mm)
      Sig1 <- mm %*% vcov1 %*% t(mm)
      # Sig0 and Sig1 are (num_sentinels x num_sentinels)

      # Standard errors
      se0 <- sqrt(var0)
      se1 <- sqrt(var1)
      se_tau <- sqrt(var_tau)

      #Add to Ys data frame:
      Ys <- data.frame(
        sentNum = sentinels$sentNum,
        rating1 = sentinels$rating1,
        rating2 = sentinels$rating2,
        Yhat0 = Yhat0,
        Yhat1 = Yhat1,
        estimate = tau_hat,
        se0 = se0,
        se1 = se1,
        se = sqrt(se1^2 + se0^2),
        se_hat = se_tau,
        weight = calc_weights( sentinels$rating1, sentinels$rating2, dat ),
        pt = paste0(round(sentinels$rating1, digits = 1), ",", round(sentinels$rating2, digits = 1))
      )

      # Ad hoc precision weights
      Ys$weight_p = 1 / Ys$se^2
      Ys$weight_p = Ys$weight_p / sum(Ys$weight_p, na.rm=TRUE)

      out_paraL <- Ys %>%
        dplyr::summarize( AFE_wt = stats::weighted.mean(estimate, w=weight, na.rm=TRUE),
                          SE_wt = sqrt(as.numeric(t(weight) %*% Sig0 %*% weight) +
                                         as.numeric(t(weight) %*% Sig1 %*% weight)),
                          AFE_prec = stats::weighted.mean( estimate, w = weight_p, na.rm=TRUE ),
                          SE_prec = sqrt(as.numeric(t(weight_p) %*% Sig0 %*% weight_p) +
                                           as.numeric(t(weight_p) %*% Sig1 %*% weight_p)),
                          n_sent_used = dplyr::n(),
                          sampsize = NA,
                          n_sentinel = mean(n_sentinel))

      out_paraL <- tibble::tibble(
        parameter = c("AFE_wt",  "AFE_prec"),
        estimate  = c(out_paraL$AFE_wt, out_paraL$AFE_prec),
        se = c(out_paraL$SE_wt,  out_paraL$SE_prec),
        n_sent_used = out_paraL$n_sent_used,
        sampsize    = out_paraL$sampsize,
        n_sentinel  = out_paraL$n_sentinel
      )

    }

    ## Parametric surface response: quad
    out_paraQ = tibble::tibble()
    if (include_PARA_QUAD) {

      sentinels <- make_sentinels( dat, n_sentinel = n_sentinel )
      df <- data.frame(rating1=sentinels$rating1, rating2=sentinels$rating2)

      #Build the model matrix (for predictions, ensures proper variable ordering):
      mm <- stats::model.matrix(~ rating1 + rating2 + I(rating1^2) + I(rating2^2) + rating1:rating2, data=df)

      dat_0 <- dplyr::filter(dat, T==0)
      dat_1 <- dplyr::filter(dat, T==1)

      M0 = estimatr::lm_robust(Y ~ rating1 + rating2 + I(rating1^2) + I(rating2^2) + rating1:rating2, data = dat_0)
      M1 = estimatr::lm_robust(Y ~ rating1 + rating2 + I(rating1^2) + I(rating2^2) + rating1:rating2, data = dat_1)

      #Get estimated coefficients:
      b0 <- stats::coef(M0)
      b1 <- stats::coef(M1)

      #Get predictions for each sentinel point:
      Yhat0 <- as.vector(mm %*% b0)  # Control model
      Yhat1 <- as.vector(mm %*% b1)  # Treated model
      tau_hat <- Yhat1 - Yhat0       # RD effect at each sentinel

      #Extract coefficient covariance matrices:
      vcov0 <- stats::vcov(M0)
      vcov1 <- stats::vcov(M1)

      #Vector containing the prediction variance at each sentinel point.
      var0 <- apply(mm, 1, function(x) as.numeric(t(x) %*% vcov0 %*% x))
      var1 <- apply(mm, 1, function(x) as.numeric(t(x) %*% vcov1 %*% x))
      # Variance of treatment effect estimate at each sentinel
      var_tau <- var0 + var1

      #IThe full covariance matrix of the predictions at all sentinels:
      Sig0 <- mm %*% vcov0 %*% t(mm)
      Sig1 <- mm %*% vcov1 %*% t(mm)
      # Sig0 and Sig1 are (num_sentinels x num_sentinels)

      # Standard errors
      se0 <- sqrt(var0)
      se1 <- sqrt(var1)
      se_tau <- sqrt(var_tau)

      #Add to Ys data frame:
      Ys <- data.frame(
        sentNum = sentinels$sentNum,
        rating1 = sentinels$rating1,
        rating2 = sentinels$rating2,
        Yhat0 = Yhat0,
        Yhat1 = Yhat1,
        estimate = tau_hat,
        se0 = se0,
        se1 = se1,
        se = sqrt(se1^2 + se0^2),
        se_hat = se_tau,
        weight = calc_weights( sentinels$rating1, sentinels$rating2, dat ),
        pt = paste0(round(sentinels$rating1, digits = 1), ",", round(sentinels$rating2, digits = 1))
      )

      # Ad hoc precision weights
      Ys$weight_p = 1 / Ys$se^2
      Ys$weight_p = Ys$weight_p / sum(Ys$weight_p, na.rm=TRUE)

      out_paraQ <- Ys %>%
        dplyr::summarize( AFE_wt = stats::weighted.mean(estimate, w=weight, na.rm=TRUE),
                          SE_wt = sqrt(as.numeric(t(weight) %*% Sig0 %*% weight) +
                                         as.numeric(t(weight) %*% Sig1 %*% weight)),
                          AFE_prec = stats::weighted.mean( estimate, w = weight_p, na.rm=TRUE ),
                          SE_prec = sqrt(as.numeric(t(weight_p) %*% Sig0 %*% weight_p) +
                                           as.numeric(t(weight_p) %*% Sig1 %*% weight_p)),
                          n_sent_used = n(),
                          sampsize = NA,
                          n_sentinel = mean(n_sentinel))

      out_paraQ <- tibble::tibble(
        parameter = c("AFE_wt",  "AFE_prec"),
        estimate  = c(out_paraQ$AFE_wt, out_paraQ$AFE_prec),
        se = c(out_paraQ$SE_wt,  out_paraQ$SE_prec),
        n_sent_used = out_paraQ$n_sent_used,
        sampsize    = out_paraQ$sampsize,
        n_sentinel  = out_paraQ$n_sentinel
      )
    }


    # Combine our results
    rs <- dplyr::bind_rows("OLS" = out_ols,
                    "bindingscore" = out_bs,
                    "frontier" = out_fr,
                    "gaussianp" = out_gp,
                    "gaussianp_res" = out_gp_res,
                    "loess" = out_loess,
                    "surflinear" = out_paraL,
                    "surfquad" = out_paraQ,
                    .id = "model" )
    rs <- tibble::remove_rownames(rs)

    rs <- rs %>%
      dplyr::distinct()
    rs
}

