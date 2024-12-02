## code to prepare `fakeData` dataset goes here


library( twodrdd )
fakeData = gen_dat_sim( sim = 1, cut1 = 0.2, cut2 = 0.3, n = 700, rho = 0.80, s = 1 , seed = 847)

names(fakeData)

head( fakeData$data )
fakeData = fakeData$data

head( fakeData )
fakeData$X = rnorm( nrow(fakeData) )


usethis::use_data(fakeData, overwrite = TRUE)
