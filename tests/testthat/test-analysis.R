
test_that("the analysis functio nworks", {

  data( fakeData )

  anal<- analysis( fakeData,
                   cut1=0, cut2=0,
                   n_sentinel = 20,
                   include_OLS = TRUE,
                   include_BINDING = TRUE,
                   include_FRONTIER = TRUE,
                   include_GP = TRUE,
                   include_GP_RES = TRUE,
                   include_loess = TRUE,
                   include_PARA_LINEAR = TRUE,
                   include_PARA_QUAD = TRUE )

  expect_true( is.data.frame(anal) )
  expect_true( all( !is.na( anal$estimate ) ) )
} )
