
# This shows how the rddapp is unstable with repeated calls.
# In the following, A1 and A2 should be the same.

rm(list=ls())

library(rddapp)
library(dplyr)
library(tibble)
library(purrr)

set.seed( 444944 )


my_frontier <- function( dat, cut1= 0, cut2= 0 ) {

    dat$c1 = dat$rating1 > cut1
    dat$c2 = dat$rating2 > cut2

    tb <- table( c1= dat$c1, c2=dat$c2 )
    stopifnot( all( tb[,2] >= 10 ) )
    stopifnot( all( tb[2,] >= 10 ) )

    # do RDD on c1 and c2
    runfun_frontier<-mrd_est(data=dat, Y~rating1 + rating2 | rating1 + rating2,
                             method="front",
                             cutpoint=c(cut1, cut2),
                             t.design=c("l","l"), boot=0) #boot=100


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

    het.tau <- runfun_frontier[["front"]][["tau_MRD"]][["est"]][2,6]
    het.se <- runfun_frontier[["front"]][["tau_MRD"]][["se"]][2,6]
    #het.se <- NA
    het.n <- runfun_frontier[["front"]][["tau_MRD"]][["obs"]][["bw"]][2]
    #2 is HET model

    ##Try Het Model for RV1
    het.tau1 <- runfun_frontier[["front"]][["tau_MRD"]][["est"]][2,4]
    het.se1 <- runfun_frontier[["front"]][["tau_MRD"]][["se"]][2,4]
    het.n1 <- runfun_frontier[["front"]][["tau_MRD"]][["obs"]][["bw"]][2]

    ##Try Het Model for RV2
    het.tau2 <- runfun_frontier[["front"]][["tau_MRD"]][["est"]][2,5]
    het.se2 <- runfun_frontier[["front"]][["tau_MRD"]][["se"]][2,5]
    het.n2 <- runfun_frontier[["front"]][["tau_MRD"]][["obs"]][["bw"]][2]

    out1<-c(tau.ate,se.ate,n.ate)
    out2<-c(tau.ate1,se.ate1,n.ate1)
    out3<-c(tau.ate2,se.ate2,n.ate2)
    out4<-c(het.tau,het.se,het.n)
    out5<-c(het.tau1,het.se1,het.n1)
    out6<-c(het.tau2,het.se2,het.n2)

    both <- (rbind(out1, out2, out3, out4, out5, out6))
    parameter<-c("ATE", "ATE1", "ATE2", "ATEhet", "ATEhet1", "ATEhet2")
    out_fr<-data.frame(parameter, both) %>%
        distinct()
    colnames(out_fr)<-c("parameter","estimate","se","sampsize")


    as_tibble( out_fr )
}





# Generate datasets
dat_orig_dets <- gen_dat_sim(sim = 5, rho = 0.8, s = 1, n = 1000,
                             cut1 = 0.3, cut2 = 0.3, seed = 615)


dat_orig <- dat_orig_dets$data

if ( FALSE ) {
    cat( "Fit A1\n" )
    A1 <- analysis(
        dat = dat_orig,
        min_sample = 8, startnum = NULL, endnum = NULL, n_sentinel = 20,
        include_OLS = FALSE, include_BINDING = FALSE, include_FRONTIER = TRUE,
        include_PARA_LINEAR = FALSE, include_PARA_QUAD = FALSE,
        include_GP = FALSE, include_GP_RES = FALSE, include_APPX_GP = FALSE,
        include_loess = FALSE
    )

    cat( "Fit A2\n" )
    A2 <- analysis(
        dat = dat_orig,
        min_sample = 8, startnum = NULL, endnum = NULL, n_sentinel = 20,
        include_OLS = FALSE, include_BINDING = FALSE, include_FRONTIER = TRUE,
        include_PARA_LINEAR = FALSE, include_PARA_QUAD = FALSE,
        include_GP = FALSE, include_GP_RES = FALSE, include_APPX_GP = FALSE,
        include_loess = FALSE
    )

    cat( "Fit A3\n" )
    A3 <- analysis(
        dat = dat_orig,
        min_sample = 8, startnum = NULL, endnum = NULL, n_sentinel = 20,
        include_OLS = FALSE, include_BINDING = FALSE, include_FRONTIER = TRUE,
        include_PARA_LINEAR = FALSE, include_PARA_QUAD = FALSE,
        include_GP = FALSE, include_GP_RES = FALSE, include_APPX_GP = FALSE,
        include_loess = FALSE
    )

} else {
    cat( "Fit A1\n" )
    A1 <- my_frontier( dat = dat_orig_dets$data )
    cat( "Fit A2\n" )
    A2 <- my_frontier( dat = dat_orig_dets$data )
    cat( "Fit A3\n" )
    A3 <- my_frontier( dat = dat_orig_dets$data )
}

print( A1 )
print( A2 )
A1$delta = A1$estimate - A2$estimate
A1$delta2 = A1$estimate - A3$estimate
A1$delta23 = A2$estimate - A3$estimate

cat( "The delta columns should all be zero:\n" )
A1 %>%
    relocate( delta, delta2, delta23, .after = estimate ) %>%
    print()


