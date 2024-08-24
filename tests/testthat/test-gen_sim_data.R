

test_that("generating simulation data", {

  library(tidyverse)

  dat = gen_dat_sim(4, 0.8, 0.75, 0.2, 444343, n = 5000, s = 1)
  dat2 = gen_dat_sim(4, 0.8, 0.75, 0.2, 444343, n = 5000, s = 1)
  expect_equal( dat$parameters, dat2$parameters )
  expect_equal( dat$data, dat2$data )

  names(dat)

  expect_true( is.function( dat$tx.func ) )
  expect_equal( dat$tx.func( 0, 0 ), 0.4 )

  dat$parameters
  expect_true( is.list(dat$parameters) )

  dat = gen.2RRDD(
    n = 5000,
    s = 1,
    mu1 = 10,
    mu2 = 5,
    sigma.r1 = 1,
    cut1.quantile = 0.5, cut2.quantile = 0.95,
    param.Y0 = c( 1,2,3,4,5,6 ),
    tx.func = function( r1, r2 ) {
      0.2 + r1
    }
  )


  dat$parameters
  dat = dat$data
  expect_equal( nrow(dat), 5000 )
  head(dat)
  expect_true( all( c("rating1", "rating2", "T", "Y") %in% names(dat) ) )
  expect_true( length( unique(dat$T) ) <= 2 )
  expect_equal( length(unique(dat$sampleID)), 1 )

  summary(dat$Y)

})
