
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
         include_BINDING = TRUE,
         include_FRONTIER = TRUE,
         include_OLS = TRUE,
         include_GP = TRUE,
         include_GP_RES = TRUE,
         include_loess = TRUE
         )

res |> dplyr::select( -n ) |>
  knitr::kable( digits = 2)
```

| model         | parameter | estimate |   se | n_sent | sampsize |
|:--------------|:----------|---------:|-----:|-------:|---------:|
| OLS           | AFE       |    -0.60 | 0.46 |     NA |      700 |
| bindingscore  | ATE       |     0.42 | 0.19 |     NA |      638 |
| bindingscore  | ATE       |     0.42 | 0.19 |     NA |      638 |
| frontier      | ATE       |     0.61 | 0.25 |     NA |      682 |
| frontier      | ATE1      |    -0.05 | 0.29 |     NA |      682 |
| frontier      | ATE2      |     1.13 | 0.37 |     NA |      682 |
| gaussianp     | AFE_wt    |     0.60 | 0.21 |     20 |       NA |
| gaussianp     | AFE_prec  |     0.59 | 0.19 |     20 |       NA |
| gaussianp_res | AFE_wt    |     0.63 | 0.21 |     20 |       NA |
| gaussianp_res | AFE_prec  |     0.61 | 0.19 |     20 |       NA |
| loess         | AFE_wt    |     1.11 |   NA |     39 |       NA |
| loess         | AFE_prec  |     0.81 |   NA |     39 |       NA |

There are other functions to get more rich descriptions of the analysis.
See, e.g., `gaussianp()`.
