## <img src = "docs/figures/bendr_hex.png" align="right" width="400" height = "400"> rstapDP: Spatial Temporal Aggregated Dirichlet Process Predictors in R
<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## About

This is an R package that fits a Normal Linear Model to data with a Dirichlet Process Prior placed on a subset of the regression coefficients.
The primary target audience is researchers interested in the effect of built environment features (BEFs) on human health, 
though other applications are possible. 

## Installation

### Development Version

 Currently this package is only available via Github. In order to install the software use the following 
 lines of R code

 ```r
 if(!require(devtools)){
	install.packages("devtools")
	library(devtools)
 }

install_github("apeterson91/bendr",dependencies = TRUE)
 ```

## Contributing

 Examples and code contributions are welcome. Feel free to start/address a feature in the issue tracker and I'll be notified shortly. 

#### Code of Conduct

Please note that `rstapDP` is released with a [Contributor Code of Conduct](https://www.contributor-covenant.org/). By contributing to this project, you agree to abide by its terms.


## How to cite this package

 A citation is in progress. Check back soon.

## Acknowledgments

This work was developed with support from NIH grant R01-HL131610 (PI: Sanchez).


