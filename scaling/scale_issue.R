# Exploring scale issue


library( tidyverse )
library( twodrdd )

set.seed( 404402 )

# Helper: run analysis and gather desired outputs for a given label and dataset
analyze_and_extract <- function(label, dat, sents = NULL) {
    cat( "analyzing", label, "\n" )
    # Run analysis
    out <- analysis(
        dat = dat,
        startnum = NULL, endnum = NULL, n_sentinel = 20,
        include_OLS = FALSE, include_BINDING = TRUE, include_FRONTIER = TRUE,
        include_PARA_LINEAR = FALSE, include_PARA_QUAD = FALSE,
        include_GP = FALSE, include_GP_RES = FALSE,
        include_loess = FALSE
    )


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
        row = label,
        sd1 = sd(dat$rating1),
        sd2 = sd(dat$rating2),
        w1 = b$R1,
        w2 = b$R2,
        length1 = mean( lengthscales$length1 ),
        length2 = mean( lengthscales$length2 ),
        AFE_wt = get_model_param(out, "gaussianp_res", "AFE_wt"),
        AFE_prec = get_model_param(out, "gaussianp_res", "AFE_prec"),
        BS_AFE = get_model_param(out, "bindingscore", "AFE"),
        PF_AFE = get_model_param(out, "frontier", "AFE"),
        PF_AFE1 = get_model_param(out, "frontier", "AFE1"),
        PF_AFE2 = get_model_param(out, "frontier", "AFE2")
    )
}




# Generate scaled datasets ----
dat_orig_dets <- gen_dat_sim(sim = 5, rho = 0.8, s = 1, n = 5000,
                             cut1 = 0.3, cut2 = 0.3, seed = 615)

dat_orig_dets$parameters

dat_orig <- dat_orig_dets$data[1:700,]

# Look at the true treatment curve
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


if ( FALSE ) {
    # testing the function
    analyze_and_extract( "og", dat_orig )
}

# Fit to sequence of scaled models ----



# Create all datasets in a named list
all_dats <- list(
    oo = dat_orig,
    original = dat_orig,
    scaled0.5 = dat_orig %>% mutate(rating1 = rating1 * 0.5),
    scaled1.1 = dat_orig %>% mutate(rating1 = rating1 * 1.1),
    scaled1.5 = dat_orig %>% mutate(rating1 = rating1 * 1.5),
    scaled2 = dat_orig %>% mutate(rating1 = rating1 * 2),
    scaled10 = dat_orig %>% mutate(rating1 = rating1 * 10),
    scaled100 = dat_orig %>% mutate(rating1 = rating1 * 100)
)


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


# Calculate for each and combine
full_result_table <- map2_dfr( names(all_dats), all_dats,
                               analyze_and_extract )

# Drop first row due to weird frontier method glitch
full_result_table <- full_result_table %>%
    filter(row != "oo") %>%
    arrange( sd1 )


# Table of results ----
knitr::kable( full_result_table, digits = 2 )


# Looking at the true DGP ----

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






