#Lily An
#Compares GPR to binding score, loess, and frontier



#' code to analyze single dataset of n observations using our list of
#' estimators.
#'
#' @param radius Parameter to loess estimator.
#' @param min_sample Parameter to loess estimator.
#' @param startnum Parameter to gaussian process estimator
#' @param dat data to be analyzed
#' @param n_sentinel number of sentinels per side
#' @param include_OLS for OLS regression
#' @param include_BINDING for binding score
#' @param include_FRONTIER for pooled frontier
#' @param include_GP for GPR
#' @param include_GP_RES for residualized GPR
#' @param include_loess for loess
#' @param endnum Parameter to gaussian process estimator.
#'
#' @export
analysis <- function( dat, radius=10, min_sample=8, n_sentinel = 20,
                      startnum=50, endnum=100,
                      include_OLS = TRUE,
                      include_BINDING = FALSE,
                      include_FRONTIER = FALSE,
                      include_GP = FALSE,
                      include_GP_RES = FALSE,
                      include_loess = FALSE ) {

    Y <- NULL

    # Assumption: ratings have been centered at cutpoints (so we can use
    # cut1 = cut2 = 0).  Our DGP centered the ratings automatically, so
    # this is ok.
    cut1 = 0
    cut2 = 0

    out_ols = tibble::tibble()
    out_bs = tibble::tibble()
    out_fr = tibble::tibble()

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


    ### Binding Score y ~ x1 + x2 | c1 + c2 for a sharp MRDD with two covariates
    #Use the package's CV optimal bandwidth
    if ( include_BINDING ) {
        runfun_bs <- rddapp::mrd_est(data=dat, Y~rating1 + rating2 | rating1 + rating2,
                             method="center",
                             cutpoint=c(cut1, cut2),
                             t.design=c("l","l"))
        #Could add: se.type="vcovHC", default is "HC1"

        tau.ate <- runfun_bs[["center"]][["tau_MRD"]][["est"]][4]
        se.ate <- runfun_bs[["center"]][["tau_MRD"]][["se"]][4]
        n.ate <- runfun_bs[["center"]][["tau_MRD"]][["obs"]][4]
        #4 is the Optimal bandwidth specification
        out_bs <- c(tau.ate, se.ate, n.ate) # 3 columns - tau, se, n

        both <- rbind(out_bs, out_bs)
        parameter <- c("ATE", "ATE")
        out_bs <- data.frame(parameter, both)
        colnames(out_bs) <- c("parameter", "estimate", "se", "sampsize")
        out_bs <- out_bs %>%
            dplyr::filter(parameter=="ATE")
    }
        ## Frontier
        #Use the package's CV optimal bandwidth
    if ( include_FRONTIER ) {

        dat$c1 = dat$rating1 > cut1
        dat$c2 = dat$rating2 > cut2

        # #Either 0 or all
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

            tau.ate <- runfun_frontier1[["est"]][4]
            se.ate <- runfun_frontier1[["se"]][4]
            n.ate <- runfun_frontier1[["obs"]][4]
            #4 is the Optimal bandwidth specification
            out_fr <- c(tau.ate, se.ate, n.ate) # 3 columns - tau, se, n

            both <- rbind(out_fr, out_fr)
            parameter <- c("ATE", "ATE")
            out_fr <- data.frame(parameter, both)
            colnames(out_fr) <- c("parameter", "estimate", "se", "sampsize")
            out_fr <- out_fr %>%
                dplyr::filter(parameter=="ATE")

         } else if ( any( tb[2,] <= 10 ) ) {
            # do RDD on c2 only
            runfun_frontier2<-rddapp::rd_est(data=dat, Y~rating2,
                                    cutpoint=cut2,
                                    t.design="l")

            tau.ate <- runfun_frontier2[["est"]][4]
            se.ate <- runfun_frontier2[["se"]][4]
            n.ate <- runfun_frontier2[["obs"]][4]
            #4 is the Optimal bandwidth specification
            out_fr <- c(tau.ate, se.ate, n.ate) # 3 columns - tau, se, n

            both <- rbind(out_fr, out_fr)
            parameter <- c("ATE", "ATE")
            out_fr <- data.frame(parameter, both)
            colnames(out_fr) <- c("parameter", "estimate", "se", "sampsize")
            out_fr <- out_fr %>%
                dplyr::filter(parameter=="ATE")

        } else {
            # do RDD on c1 and c2
        runfun_frontier<-rddapp::mrd_est(data=dat, Y~rating1 + rating2 | rating1 + rating2,
                                 method="front",
                                 cutpoint=c(cut1, cut2),
                                 t.design=c("l","l"), boot=50)

        ##Non-parametric using optimal bandwidth case (bw selected by CV)
        tau.ate <- runfun_frontier[["front"]][["tau_MRD"]][["est"]][2,3]
        se.ate <- runfun_frontier[["front"]][["tau_MRD"]][["se"]][2,3]
        #se.ate <- NA
        n.ate <- runfun_frontier[["front"]][["tau_MRD"]][["obs"]][["bw"]][1]
        #1 is constant effects model model

        ##Try Complete Model for RV1
        ##Non-parametric using optimal bandwidth case (bw selected by CV)
        tau.ate1 <- runfun_frontier[["front"]][["tau_MRD"]][["est"]][2,1]
        se.ate1 <- runfun_frontier[["front"]][["tau_MRD"]][["se"]][2,1]
        n.ate1 <- runfun_frontier[["front"]][["tau_MRD"]][["obs"]][["bw"]][1]

        ##Try Complete Model for RV2
        ##Non-parametric using optimal bandwidth case (bw selected by CV)
        tau.ate2 <- runfun_frontier[["front"]][["tau_MRD"]][["est"]][2,2]
        se.ate2 <- runfun_frontier[["front"]][["tau_MRD"]][["se"]][2,2]
        n.ate2 <- runfun_frontier[["front"]][["tau_MRD"]][["obs"]][["bw"]][1]

        out1<-c(tau.ate,se.ate,n.ate)
        out2<-c(tau.ate1,se.ate1,n.ate1)
        out3<-c(tau.ate2,se.ate2,n.ate2)

        both <- (rbind(out1, out2, out3))
        parameter<-c("ATE", "ATE1", "ATE2")
        out_fr<-data.frame(parameter, both) %>%
            dplyr::distinct()
        colnames(out_fr)<-c("parameter","estimate","se","sampsize")
        }
    }

    ## GaussianP
    #safe_gp = safely(gaussianp)
    out_gp = tibble::tibble()
    if ( include_GP ) {
      runfun_gaussianp <- gaussianp(sampdat=dat,
                                    n_sentinel = n_sentinel,
                                    startnum = startnum, endnum=endnum,
                                    method = "new" )

      out_gp = calculate_average_impact( runfun_gaussianp )
      out_gp$n = nrow( dat )

    }

    out_gp_res = tibble::tibble()
    if (include_GP_RES) {
        M = stats::lm( Y ~ rating1 + rating2, data=dat )
        #Ms = summary(M)

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
        out_loess <- loess2DRDD(sampdat = dat,
                                radius = radius,
                                n_sentinel = n_sentinel )
        out_loess = calculate_average_impact( out_loess, calc_SE=FALSE )
        out_loess$n = nrow(dat)
    }

    # Combine our results
    rs <- dplyr::bind_rows("OLS" = out_ols,
                    "bindingscore" = out_bs,
                    "frontier" = out_fr,
                    "gaussianp" = out_gp,
                    "gaussianp_res" = out_gp_res,
                    "loess" = out_loess,
                    .id = "model" )
    rs <- tibble::remove_rownames(rs)

    rs
}

