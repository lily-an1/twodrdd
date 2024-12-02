## code to prepare `fakeData` dataset goes here


library( twodrdd )
cut1 = 0.2
cut2 = 0.3
fakeData = gen_dat_sim( sim = 1, cut1 = cut1, cut2 = cut2, n = 700, rho = 0.80, s = 1 , seed = 847)

names(fakeData)

head( fakeData$data )
fakeData = fakeData$data

fakeData$rating1 = fakeData$rating1-cut1
fakeData$rating2 = fakeData$rating2-cut2

usethis::use_data(fakeData, overwrite = TRUE)
