




if ( FALSE ) {
  dat = tibble( rating1 = runif( 1000, -1, 1 ),
                rating2 = runif( 1000, -1, 1 ) )

  dat$rating1[1] = 3
  dat
  sents = twodrdd:::make_sentinels( dat, 1000,
                                    add_weights = TRUE,
                                    method = "kernel" )
  #%>%
  #  drop_low_weight_sentinels()
  sents %>%
    knitr::kable( digits = 3 )

  sum(sents$weight)
  sum( sents$weight[ sents$rating1 == 0 ] )
  sum( sents$weight[ sents$rating2 == 0 ] )
}


dat = tibble( rating1 = runif( 10000, -1.5, 0.5 ),
              rating2 = runif( 10000, -1, 1 ) )



sents = twodrdd:::make_sentinels( dat, 20,
                                  add_weights = TRUE, method="kernel" )
sents



calc_edge_weight <- function( sents ) {
  tot <- sum(sents$weight)
  R2 = sum( sents$weight[ sents$rating1 == 0 ] )
  R1 = sum( sents$weight[ sents$rating2 == 0 ] )

  tibble( total = tot, R1 = R1, R2 = R2 )
}


ggplot( dat, aes( rating1, rating2 ) ) +
  geom_point() +
  geom_point( data=sents, aes( size =weight ), col="red" ) +
  coord_fixed()


ggplot( sents, aes( sentNum, weight ) ) +
  geom_point()

calc_edge_weight( sents )


dat2 = dat
dat2$rating1 = dat2$rating1 * 10
sents2 = twodrdd:::make_sentinels( dat2, 20,
                                   add_weights = TRUE, method="kernel" )
sents2
ggplot( sents2, aes( sentNum, weight ) ) +
  geom_point()
calc_edge_weight( sents2 )



dat3 = dat
dat3$rating1[1] = 3
dat3
sents3 = twodrdd:::make_sentinels( dat3, 20,
                                   add_weights = TRUE, method="kernel" )
ggplot( sents3, aes( sentNum, weight ) ) +
  geom_point()

sents3 %>%
  print( n = 100 )

drop_low_weight_sentinels( sents3 ) %>%
  knitr::kable( digits = 3 )

diff(sents3$rating1 )
diff(sents3$rating2 )

ggplot( dat3, aes( rating1, rating2 ) ) +
  geom_point( alpha=0.1 ) +
  geom_point( data=sents3, aes( size =weight ), col="red" ) +
  coord_fixed()

calc_edge_weight( sents3 )



