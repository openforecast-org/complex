# complex package for R
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/complex)](https://cran.r-project.org/package=complex)
[![Downloads](https://cranlogs.r-pkg.org/badges/complex)](https://cran.r-project.org/package=complex)
[![R-CMD-check](https://github.com/openforecast-org/complex/actions/workflows/test.yml/badge.svg)](https://github.com/openforecast-org/complex/actions/workflows/test.yml)

Time series analysis and forecasting using complex variables

![hex-sticker of the complex package for R](https://github.com/openforecast-org/complex/blob/master/man/figures/complex-web.png?raw=true)

The package includes basic instruments for correlation and regression analysis of complex-valued variables. The package supports the monograph by Svetunkov & Svetunkov "Complex-valued Econometrics with Examples in R" <[doi:10.1007/978-3-031-62608-1](https://doi.org/10.1007%2F978-3-031-62608-1)>.


### Installation
Stable version can be installed from CRAN:
```r
install.packages("complex")
```

For installation from github use the `remotes` package:
```r
if (!require("remotes")){install.packages("remotes")}
remotes::install_github("openforecast-org/complex")
```
