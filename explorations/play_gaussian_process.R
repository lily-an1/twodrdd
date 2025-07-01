

library( tidyverse )
library(MASS)

set.seed( 404404 )

simulate_gp <- function(n = 700, lengthscales = c(0.5, 1), nugget = 1) {
  X <- matrix(runif(n * 2, -2, 2), ncol = 2)

  # Anisotropic squared exponential kernel
  sq_exp_kernel <- function(x1, x2, ls) {
    dists <- sweep(x1, 2, ls, "/") - sweep(x2, 2, ls, "/")
    d2 <- rowSums(dists^2)
    exp(-0.5 * d2)
  }

  K <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] <- exp(-sum((X[i, ] - X[j, ])^2 / (2 * lengthscales^2)))
    }
  }

  diag(K) <- diag(K) + nugget

  Y <- mvrnorm(1, mu = rep(0, n), Sigma = K)

  tibble(x1 = X[,1], x2 = X[,2], y = Y)
}



library(laGP)

# Current thoughts: theta (lengthscales) need to be square rooted to
# be put back in linear terms.
#
# The nugget seems to be something like signal to noise, and is
# implicitly scaling Y to unit variance.  That said, rescaling Y does
# change the value of the estimated nugget, so something weird is
# going on in the fitting procedure.

# Simulate data
dat <- simulate_gp()
#dat$y = dat$y / sd( dat$y )#* 100
dat$x1 = dat$x1


# Plot the data
library(akima)

interp_dat <- with(dat, akima::interp(x1, x2, y, duplicate = "mean"))

interp_df <- expand.grid(x = interp_dat$x, y = interp_dat$y)
interp_df$z <- as.vector(interp_dat$z)

ggplot(interp_df, aes(x = x, y = y, fill = z)) +
  geom_raster() +
  geom_point(data = dat, aes(x = x1, y = x2), inherit.aes = FALSE,
             color = "black", size = 0.5) +
  scale_fill_viridis_c() +
  theme_minimal()



X <- as.matrix(dat[, c("x1", "x2")])
Y <- dat$y

# Set up priors
da <- darg(list(mle = TRUE, isotropic = FALSE), X)
ga <- garg(list(mle = TRUE), Y)

# Fit GP
gpi <- newGPsep(X, Y, d = da$start, g = ga$start, dK = TRUE)

mle <- mleGPsep(gpi, param = "both", tmin = c(da$min, ga$min), tmax = c(da$max, ga$max),
         ab = c(da$ab, ga$ab), maxit = 1000)

print( mle )

