
test_that("sentinel code works", {


  data( fakeData )

  fakeData$rating1 = runif( nrow(fakeData), min = -1, max = 1 )

  sents <- twodrdd:::make_sentinels(fakeData, n_sentinel = 10, fixed_sent = FALSE)
  sents

  sents$w1 = twodrdd:::calc_weights( sents$rating1, sents$rating2, data=fakeData )
  sents$w2 = twodrdd:::calc_weights( sents$rating1, sents$rating2, data=fakeData, method = "kernel" )

  expect_true( all( sents$w1 >= 0 ) )
  expect_true( all( sents$w2 >= 0 ) )

  if ( FALSE ) {
    sL <- pivot_longer(sents, cols = starts_with("w"),
                        names_to = "weight_type", values_to = "weight" )
    ggplot( sL, aes( sentNum, weight, col=weight_type ) ) +
      geom_point() + geom_line()


    ggplot( fakeData, aes(rating1, rating2, col=T ) ) +
      geom_point()

    dat2 = fakeData
    dat2$rating1 = dat2$rating1 * 10
    sents2 <- twodrdd:::make_sentinels(dat2, n_sentinel = 10, fixed_sent = FALSE)
    sents2$w1 = twodrdd:::calc_weights( sents2$rating1, sents2$rating2, data=dat2 )
    sents2$w2 = twodrdd:::calc_weights( sents2$rating1, sents2$rating2, data=dat2, method = "kernel" )

    sL2 <- pivot_longer(sents2, cols = starts_with("w"),
                        names_to = "weight_type", values_to = "weight" )
    ggplot( sL2, aes( sentNum, weight, col=weight_type ) ) +
      geom_point() + geom_line()

    # weights are unchanged by scaling
    sL2$weight / sL$weight
  }
})
