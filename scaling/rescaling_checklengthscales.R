#figure 8 in paper

# Is rescaling a problem for 2D RDDs? Is it unique to 2D RDDs?

## DGP (assuming two lengthscales):
## No TE on one boundary, growing TE on other boundary,

library( twodrdd )

library( tidyverse )



if ( FALSE ) {

  dat = gen_dat_sim(sim=5, rho=0.8, s = 1, n = 700)

  ggplot( dat$data, aes( rating1, Y ) ) +
    geom_point()

  dat_orig <- gaussianp(sampdat = dat$data,
                        method="new")
  dat_orig
  params( dat_orig )

  dat2 = dat$data
  dat2$Y[dat2$T==1] = dat2$Y[dat2$T==1] + rnorm( nrow(dat2[dat2$T==1,]), sd=10 )

  dat_orig2 <- gaussianp(sampdat = dat2,
                        method="new")
  dat_orig2
  params( dat_orig2 )
  params( dat_orig )
}



#change "method" to iso if want 1 lengthscale, new for 2 lengthscales
one_run_LS <- function(sim=5, rho=0.8, s=1, n=1000, scale = 100 ) {

  dat = gen_dat_sim(sim=sim, rho=rho, s = s, n = n)
  dat_orig <- gaussianp(sampdat = dat$data,
                           method="new")

  lengthscales = params( dat_orig )


  # Rescale the data set
  newdat <- dat$data %>%
    mutate(rating1 = rating1*scale)

  # Then refit it, and show differences in length scales
  dat_rescaled <- gaussianp(sampdat = newdat,
                            method="new")
  lengthscales2 = params( dat_rescaled )

  all_length = bind_rows( orig = lengthscales, scaled = lengthscales2,
                          .id = "scaling" ) %>%
    dplyr::select( scaling, group:nugget ) %>%
    pivot_wider( names_from=group,
                 values_from=length1:nugget )

  all_length


  avgImp = calculate_average_impact(dat_orig, calc_SE = TRUE)
  avgImp_rescaled = calculate_average_impact(dat_rescaled, calc_SE = TRUE)

  all_length$AFE_wt = c(avgImp$estimate[[1]],
                        avgImp_rescaled$estimate[[1]])
  all_length$se_wt = c(avgImp$se[[1]],
                        avgImp_rescaled$se[[1]])
  all_length$AFE_prec = c(avgImp_rescaled$estimate[[2]],
                          avgImp_rescaled$estimate[[2]])
  all_length$se_prec = c(avgImp$se[[2]],
                        avgImp_rescaled$se[[2]])

  return(all_length)

}

if ( FALSE ) {
  test = one_run_LS(sim=5, rho=0.8, s=1, n=700)
}


# Run the function 10 times and add run IDs
# FLAG - I think this took a little bit
results_df <- map_dfr(c( 0.5, 1.5, 2, 5, 10, 50, 100),
                      function(scale) {
  result <- one_run_LS( scale = scale )  # Call your function
  mutate(result, scale = scale )  # Add the run ID
})

results_df <- results_df %>%
  mutate(rescaled = rep(c("original data","rescaled data"), times=10))

#Analzye
results_df %>%
  group_by(rescaled) %>%
  summarize(mean_origC_1 = mean(LS0_1),
            mean_origC_2 = mean(LS0_2),
            mean_origT_1 = mean(LS1_1),
            mean_origT_2 = mean(LS1_2),
            sd_origC_1 = sd(LS0_1),
            sd_origC_2 = sd(LS0_2),
            sd_origT_1 = sd(LS1_1),
            sd_origT_2 = sd(LS1_2))
# save(results_df, file="C:/Users/lilya/OneDrive - Harvard University/HGSE PhD/CARES Lab/2D RDD/JEBS submission/revisions/rev_lengthscale.RData")
# load("C:/Users/lilya/OneDrive - Harvard University/HGSE PhD/CARES Lab/2D RDD/JEBS submission/revisions/rev_lengthscale.RData")

windowsFonts(Times=windowsFont("Times New Roman"))
# control side
control <- ggplot(data = results_df, aes(LS0_1, LS0_2)) +
  geom_point() +
  facet_wrap(vars(rescaled)) +
  labs(x = "Lengthscale 1",
       y= "Lengthscale 2",
       subtitle = "Control GP") +
  theme(text=element_text(family="Times"))
# treated side
treat <- ggplot(data = results_df, aes(LS1_1, LS1_2)) +
  geom_point() +
  facet_wrap(vars(rescaled)) +
  labs(x = "Lengthscale 1",
       y= "Lengthscale 2",
       subtitle = "Treated GP") +
  theme(text=element_text(family="Times"))

library(patchwork)
control/treat
