# Exploring scale issue


library( tidyverse )
library( twodrdd )

set.seed( 404402 )


# Helper: run analysis and gather desired outputs for a given label and dataset
analyze_and_extract <- function(dat, scaling, sents = NULL, include_GP = TRUE ) {
    cat( "analyzing scale ", scaling , "\n" )

    dat <- dat %>%
        mutate(rating1 = rating1 * scaling)

    # Run analysis
    out <- analysis(
        dat = dat,
        startnum = NULL, endnum = NULL, n_sentinel = 20,
        include_OLS = FALSE, include_BINDING = TRUE, include_FRONTIER = TRUE,
        include_PARA_LINEAR = FALSE, include_PARA_QUAD = FALSE,
        include_GP = FALSE, include_GP_RES = FALSE,
        include_loess = FALSE,
        boot = NULL
    )


    if ( include_GP ) {
        outgp <- gaussianp( sampdat = dat,
                            startnum = NULL, endnum = NULL, n_sentinel = 20,
                            method="new", residualize = TRUE )

        out_gp_res = calculate_average_impact( outgp )
        out_gp_res$n = nrow( dat )

        lengthscales = params( outgp )
        out_gp_res$lengthscale1 = mean( sqrt( lengthscales$length1 ) )
        out_gp_res$lengthscale2 = mean( sqrt( lengthscales$length2 ) )
        out_gp_res$model = "GP_res"
        out = bind_rows( out, out_gp_res )
    } else {
        lengthscales = list( length1 = NA, length2 = NA )
    }

    # Recreate the sentinels for assessment
    if (is.null(sents)) {
        sents <- make_sentinels(dat, n_sentinel = 20,
                                fixed_sent = FALSE,
                                add_weights = TRUE,
                                drop_low_weight = TRUE )
    }

    # Boundary weights
    b = calc_edge_weight(sents)

    # Helper extractors
    get_model_param <- function(out, modelname, param) {
        value <- out %>%
            filter(model == modelname, parameter == param) %>%
            pull(estimate)
        if (length(value) == 0) NA else value[[1]]
    }

    # Assemble results
    tibble(
        scaling = scaling,
        sd1 = sd(dat$rating1),
        sd2 = sd(dat$rating2),
        w1 = b$R1,
        w2 = b$R2,
        length1 = mean( lengthscales$length1 ),
        length2 = mean( lengthscales$length2 ),
        GP_wt = get_model_param(out, "GP_res", "AFE_wt"),
        GP_prec = get_model_param(out, "GP_res", "AFE_prec"),
        BS_AFE = get_model_param(out, "bindingscore", "AFE"),
        PF_AFE = get_model_param(out, "frontier", "AFE"),
        PF_AFE1 = get_model_param(out, "frontier", "AFE1"),
        PF_AFE2 = get_model_param(out, "frontier", "AFE2")
    )
}




# Generate scaled datasets ----
dat_orig_dets <- gen_dat_sim(sim = 5, rho = 0.8, s = 1, n = 5000,
                             cut1 = 0.5, cut2 = 0.5, seed = 615)

dat_orig_dets$parameters$effects

dat_orig <- dat_orig_dets$data[1:700,]



# Look at the true treatment curve and data ----

if ( FALSE ) {

    sents <- make_sentinels(dat_orig, n_sentinel = 20, fixed_sent = FALSE,
                            add_weights = TRUE,
                            drop_low_weight = TRUE )


    tx.func <- dat_orig_dets$tx.func
    sents$tau = dat_orig_dets$tx.func( sents$rating1, sents$rating2 )
    ggplot( sents, aes( sentNum, tau ) ) +
        geom_line( aes( y = tx.func(rating1, rating2) ), col="red" ) +
        geom_point( aes( size = weight ) ) +
        labs(title = "Sentinel Estimates vs True Treatment Effect",
             x = "Sentinel Number", y = "Estimated Treatment Effect") +
        scale_size_continuous(name = "Weight") +
        theme_minimal() +
        geom_hline(yintercept = 0, linetype = "dashed", color = "blue")


    dat_orig_dets <- gen_dat_sim(sim = 5, rho = 0.8, s = 1, n = 5000,
                                 cut1 = 0.5, cut2 = 0.5, seed = 615)

    dat_orig_dets$parameters$effects
    dd = dat_orig_dets$data

    sd( dd$Y )
    summary( dd$Y )

    ggplot( dd, aes( rating1, rating2, col=as.factor(T) ) ) +
        geom_point( alpha=0.2) +
        labs(title = "Scatter of Ratings 1 and 2 by T",
             x = "Rating 1", y = "Rating 2") +
        theme_minimal()

    #dd$Y[dd$T==1] = dd$Y[dd$T==1] + 40
    ggplot( filter( dd, rating2 > 0 ), aes( rating1, Y, col=as.factor(T) ) ) +
        geom_point( alpha=0.2) +
        geom_smooth( aes( group = T ), se = FALSE, col="black" )

    ggplot( filter( dd, rating1 > 0 ), aes( rating2, Y, col=as.factor(T) ) ) +
        geom_point( alpha=0.2) +
        geom_smooth( aes( group = T ), se = FALSE, col="black" )


}


if ( FALSE ) {
    # testing the function
    sdd = as.numeric( Sys.time() )
    sdd
    dat_orig <- gen_dat_sim(sim = 5, rho = 0.8, s = 1, n = 1000,
                            cut1 = 0.5, cut2 = 0.5, seed = sdd )

    one_res <- analyze_and_extract( dat_orig$data, scaling = 1,
                                    include_GP = FALSE )
    one_res
}


# Fit to sequence of scaled models ----

scalings = c( 1, 1, 0.5, 1.1, 1.5, 2, 10, 100 )



one_sequence <- function( seed ) {

    dat_orig_dets <- gen_dat_sim(sim = 5, rho = 0.8, s = 1, n = 700,
                                 cut1 = 0.5, cut2 = 0.5, seed = seed)

    res = map_df( scalings,
                  \( scaling ) {
                      analyze_and_extract( dat_orig_dets$data, scaling = scaling )
                  } )

    res$seed = seed

    res[ -1, ]
}


if ( FALSE ) {
    one_sequence( 343434 )
}

cli::cli_h1("Starting the simulation" )

library( furrr )
library( future )
nwork = parallel::detectCores() - 2
nwork
plan(multisession, workers = nwork )

R = ceiling( 100 / nwork )
R

tictoc::tic()
res = future_map( 17 * 1:(R*nwork) + 2235523,
                  .f = one_sequence,
                  .options = furrr_options(seed = NULL),
                  .progress = TRUE )
print( tictoc::toc() )
plan(sequential)

cli::cli_h1("Combining results" )
res <- bind_rows( res, .id ="runID" )
res

saveRDS( res, here::here( "scaling/scaling_results.rds" ) )

res = read_rds( here::here( "scaling/scaling_results.rds" ) )

# Look at variation in estimates
sums <- res %>% group_by( scaling ) %>%
    summarise( across( c(w1, w2, length1, length2,
                         GP_wt, GP_prec, BS_AFE, PF_AFE, PF_AFE1, PF_AFE2),
                       c( mn = mean, sd = sd ) ) )
sums

# Make table for the supplement of the paper
sums2 <- res %>% group_by( scaling ) %>%
    summarise( across( c(w1, w2, length1, length2,
                         GP_wt, GP_prec, BS_AFE, PF_AFE, PF_AFE1, PF_AFE2),
                       mean ) )
sums2

sums2 %>%
    dplyr::select( -PF_AFE1, -PF_AFE2 ) %>%
xtable::xtable( ) %>%
    print( include.rownames = FALSE )



# Look at how weights change with scaling

if ( FALSE ) {
    asents <- map( all_dats, function( ddd ) {
        make_sentinels(ddd, n_sentinel = 20, fixed_sent = FALSE) %>%
            mutate(weight = calc_weights(rating1, rating2, ddd)) %>%
            drop_low_weight_sentinels()
    } )

    asents[[1]]
    asents[[4]]


    sent1 = make_sentinels( all_dats$original, n_sentinel = 20, fixed_sent = FALSE)
    sent100 = make_sentinels( all_dats$scaled100, n_sentinel = 20, fixed_sent = FALSE)
    tail( sent1 )
    tail( sent100 )

    calc_weights( sent1$rating1, sent1$rating2, all_dats$original )


    analyze_and_extract( "tst", all_dats$original )
    analyze_and_extract( "tst", all_dats$scaled1.1 )

}


# Looking at the true DGP ----

if ( FALSE ) {
    dat_orig_dets$parameters$effects


    d_big <- dat_orig_dets$data
    F1 = d_big %>%
        filter( abs( rating1 ) < 0.25, rating2 > 0 )
    F1 %>%
        group_by( T ) %>%
        summarise( meanY = mean(Y, na.rm = TRUE),
                   n = n(),
                   SE = sd(Y)/sqrt( n() ) )

    ggplot( F1, aes( rating2, Y, col = as.factor(T) ) ) +
        geom_point() +
        geom_smooth( method = "loess", se = FALSE ) +
        labs(title = "Mean Y by T for rating1 near zero",
             x = "Rating 2", y = "Y") +
        theme_minimal()


    F2 = d_big %>%
        filter( abs( rating2 ) < 0.25, rating1 > 0 )
    F2 %>%
        group_by( T ) %>%
        summarise( meanY = mean(Y, na.rm = TRUE),
                   n = n(),
                   SE = sd(Y)/sqrt( n() ) )

    ggplot( F2, aes( rating1, Y, col = as.factor(T) ) ) +
        geom_point() +
        geom_smooth( method = "loess", se = FALSE ) +
        labs(title = "Mean Y by T for rating1 near zero",
             x = "Rating 2", y = "Y") +
        theme_minimal()

    ggplot( dat_orig, aes( rating1, rating2, col = as.factor(T) ) ) +
        geom_point() +
        labs(title = "Scatter of Ratings 1 and 2 by T",
             x = "Rating 1", y = "Rating 2") +
        theme_minimal()

    dat_orig_dets$parameters$effects
}






