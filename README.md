
# Package to fit 2D RDD models

This package fits 2D RDD models using loess or gaussian processes. The
package is still in development and is not yet on CRAN.

## Installation

You can install the development version of twodrdd from
[GitHub](https://github.com/lily-an1/twodrdd) with:

    # install.packages("devtools")
    devtools::install_github("lily-an1/twodrdd")

## Example

This is a basic example which shows you how to fit a 2D RDD model:

``` r
library( twodrdd )
dat = gen_dat_sim( sim = 5, n = 700, rho = 0.80, s = 1 )
res <- analysis(dat$data,
         n_sentinel = 20,
         include_BINDING = TRUE,
         include_FRONTIER = TRUE,
         include_OLS = TRUE,
         include_GP = TRUE,
         include_GP_RES = TRUE,
         include_loess = TRUE,
         include_PARA_LINEAR = FALSE,
         include_PARA_QUAD = FALSE
         )

knitr::kable( res, digits = 2)
```

| model         | parameter | estimate |   se | n_sent |   n | sampsize |
|:--------------|:----------|---------:|-----:|-------:|----:|---------:|
| OLS           | AFE       |     3.60 | 0.54 |     NA | 700 |      700 |
| bindingscore  | AFE       |     0.69 | 0.21 |     NA |  NA |      692 |
| bindingscore  | AFE       |     0.69 | 0.21 |     NA |  NA |      692 |
| frontier      | AFE       |     0.48 | 0.19 |     NA |  NA |      697 |
| frontier      | AFE1      |     0.42 | 0.28 |     NA |  NA |      697 |
| frontier      | AFE2      |     0.57 | 0.24 |     NA |  NA |      697 |
| gaussianp     | AFE_wt    |     0.36 | 0.21 |     18 | 700 |       NA |
| gaussianp     | AFE_prec  |     0.44 | 0.19 |     18 | 700 |       NA |
| gaussianp_res | AFE_wt    |     0.36 | 0.21 |     18 | 700 |       NA |
| gaussianp_res | AFE_prec  |     0.44 | 0.19 |     18 | 700 |       NA |
| loess         | AFE_wt    |     0.29 |   NA |     39 | 700 |       NA |
| loess         | AFE_prec  |     0.50 |   NA |     39 | 700 |       NA |

There are other functions to get more rich descriptions of the analysis.
See, e.g., `gaussianp()`.
