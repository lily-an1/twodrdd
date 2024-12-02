

#' \code{twodrdd} package
#'
#'
#' @references
#' An et. al (2024) Our paper on using Gaussian process regression as a nonparametric surface response method for two-dimensional RDDs in education.
#'
#' @name twodrdd
#' @importFrom dplyr %>%
#' @importFrom purrr %||%
"_PACKAGE"



#'@title Fake 2D RDD dataset
#'
#'@description This is a fake dataset for illustrative purposes for the twodrdd package.
#'
#'@format A data frame with 1500 rows and 7 variables: \describe{
#'  \item{\code{sampleID}}{integer ID number of which run is used. This is helpful in the event of conducting multiple runs in a simulation.}
#'  \item{\code{rating1}}{double Value of the first running variable, centered around its mean.}
#'  \item{\code{rating2}}{double Value of the second running variable, centered around its mean.}
#'  \item{\code{Y}}{double Observed outcome (think test score).)}
#'  \item{\code{T}}{integer Treatment assignment (1 treated, 0 control).} }
#'@details These data were generated using the `gen_dat_sim()` package.
"fakeData"
