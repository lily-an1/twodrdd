


test_that( "support functions of gaussianp work", {

  dat = gen_dat_sim( sim = 2, n = 1700, rho = 0.80, s = 1 )
  expect_true( nrow(dat$data) == 1700 )

  sampdat = dat$data
  expect_true( is.data.frame(sampdat) )

  #Take the min. rating value (score) from each row
  r.c<-apply(sampdat[,c("rating1","rating2")],1,min)

  # If the min. is below 0, mark as treated
  # note: this assumes RVs are centered around their respective cut points.
  T.c<-1*(r.c<0)

  # create dataset with outcome, "treated" (where treated if score<0),
  # what the score was, rating 1, rating 2
  sampdat<-data.frame( Y=sampdat$Y,
                       T.c,
                       r.c,
                       rating1=sampdat$rating1,
                       rating2=sampdat$rating2 )

  n = nrow(sampdat)
  rng = range(sampdat$Y )
  rng

  expect_true( -100 <= rng[[1]] & rng[[2]] <= 100 )

  sentinels = twodrdd:::make_sentinels(sampdat, n_sentinel = 10, fixed_sent = TRUE )
  expect_true( is.data.frame( sentinels ) )
  expect_equal( nrow(sentinels), 2*10 - 1 )

  sentinels <- sentinels %>%
    dplyr::mutate(weight = twodrdd:::calc_weights( sentinels$rating1,
                                         sentinels$rating2,
                                         data = sampdat ))

  expect_true( all( !is.na( sentinels$weight ) ) )

  s_dp <- sentinels %>%
    drop_low_weight_sentinels( drop_prop = 0.2 )
  expect_equal( colnames(s_dp), colnames(sentinels) )
  expect_true( nrow(s_dp) < nrow(sentinels) )
  t_wt = sum( s_dp$weight )
  expect_true( t_wt > 1 - 0.2 )
  expect_true( t_wt < 1 - 0.2/2 )

})




test_that("gaussianp works", {


  dat = gen_dat_sim( sim = 5, n = 700, rho = 0.80, s = 1 )
  head( dat$data )
  if ( FALSE ) {
    library(tidyverse)
    ggplot( dat$data, aes(rating1, rating2, col=T ) ) +
      geom_point()
  }

  dat$data$Y[ dat$data$T== 1] = dat$data$Y[ dat$data$T== 1] + 2



  # takes time to run, but running it on approximate to be faster.
  rsG_full = gaussianp(dat$data, n_sentinel = 40)
  expect_true( all( c("sentNum", "rating1", "rating2", "Yhat0", "estimate") %in% colnames(rsG_full) ) )
  expect_true( all( !is.na( rsG_full$var1.37 ) ) )

  rsG_full

  afe <- calculate_average_impact( rsG_full, calc_SE = TRUE )
  afe
  expect_true( all( afe$estimate > 0 ) )

  dat$data$Y[ dat$data$T== 1] = dat$data$Y[ dat$data$T== 1] - 3
  rsG_full = gaussianp(dat$data, n_sentinel = 20)
  expect_true( all( !is.na( rsG_full$var1.13 ) ) )

  afe <- calculate_average_impact( rsG_full, calc_SE = TRUE )
  afe
  expect_true( all( afe$estimate < 0 ) )



  if (FALSE ) {
    nrow(rsG_full)
    rsG <- drop_low_weight_sentinels( rsG_full )

    expect_true( nrow( rsG ) == nrow( rsG_full ) )

    # Should this be fewer
    expect_true( ncol( rsG ) == ncol( rsG_full ) )
    expect_true( all( !is.na( rsG$var1.37 ) ) )

    AFE <- dat$parameters$effects["AFE"]
    AFE

    rs$true = dat$tx.func( rs$rating1, rs$rating2 )

    head( rs )
    ggplot(rs, aes(sentNum, true)) +
      geom_line(aes(col = "tau")) +
      geom_line(aes(y = estimate, col = "Loess")) +
      geom_line(data=rsG, aes(y = estimate, col = "GP" ) )
    scale_x_continuous(breaks = rs$sentNum, labels = rs$pt) +
      coord_cartesian( ylim = c(-2,5) ) +
      theme(axis.text.x = element_text(angle = 90))
  }
})
