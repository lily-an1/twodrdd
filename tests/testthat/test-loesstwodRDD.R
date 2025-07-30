

test_that("loess estimator works", {

  set.seed(12342)
  dat = gen_dat_sim( sim = 5, rho = 0.80, s = 1 )

  radius = sd( dat$data$rating1 ) * 0.3
  radius

  sentinels <- make_sentinels( dat$data, n_sentinel = 40 )

  sentinels <- sentinels %>%
    dplyr::mutate(weight = calc_weights( sentinels$rating1,
                                         sentinels$rating2,
                                         data = dat$data )) %>%
    drop_low_weight_sentinels()

  rs <- loess2DRDD( dat$data, radius = radius, n_sentinel = 40,
                    sentinels = sentinels)

  expect_true( is.data.frame( rs ) )

  head( rs )

  avg <- calculate_average_impact( rs, calc_SE = FALSE )
  expect_true( is.data.frame( avg ) )
  expect_true( all( avg$estimate >= dat$parameters$effects["AFE"] - 0.2 ) )
  expect_true( all( avg$estimate <= dat$parameters$effects["AFE"] + 0.2 ) )



  # Try different bandwidths to see how the lines are shakey or not.
  # This will also make sure we don't crash under different radius
  # selections.
  Y_r1 = loess2DRDD( dat$data, radius = 1, n_sentinel = 40 ,
                     sentinels = sentinels)
  Y_r1$tau.hat.1 = Y_r1$estimate

  Y_r3 = loess2DRDD(dat$data, radius = 3, n_sentinel = 40,
                    sentinels = sentinels)
  Y_r1$tau.hat.3 = Y_r3$estimate

  Y_r05 = loess2DRDD(dat$data, radius = 0.5, n_sentinel = 40,
                     sentinels = sentinels)
  Y_r1$tau.hat.05 = Y_r05$estimate

  Y_r015 = loess2DRDD(dat$data, radius = 0.15, n_sentinel = 40,
                      sentinels = sentinels)
  Y_r1$tau.hat.015 = Y_r015$estimate

  Y_r1$true = dat$tx.func(  rs$rating1, rs$rating2 )

  head( Y_r1 )

  expect_true( is.data.frame( Y_r1 ) )

  # Plot the different curves
  if ( FALSE ) {
    Y_r1 %>% pivot_longer(cols = starts_with("tau"),
                          names_to = "quantity",
                          values_to = "value") %>%
      ggplot(aes(sentNum, value, col = quantity)) +
      geom_line( aes(y=true,col="tau"), lwd=2, alpha=0.5 ) +
      geom_line() +
      scale_x_continuous(breaks = Y_r1$sentNum, labels = Y_r1$pt) +
      coord_cartesian( ylim = c(-2,5) ) +
      theme(axis.text.x = element_text(angle = 90))
  }


})
