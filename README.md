## <img src="doc/figures/rstapDP_hex.png" align="right" width=95 height=100 /> `rstapDP`: Spatial Temporal Aggregated Dirichlet Process Predictors in R
<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Travis build status](https://travis-ci.org/apeterson91/rstapDP.svg?branch=master)](https://travis-ci.org/apeterson91/rstapDP)
<!-- badges: end -->

## About

The rstapDP package offers functions that fit a Linear Model with  a Dirichlet Process Prior 
placed on a subset of the regression coefficients used to model the nonlinear heterogeneous effect of built environment feature (BEF) exposure across space/time. 
The primary target audience is researchers interested in the effect of BEFs on human health, though other applications are possible. 

## Installation

### Development Version

 Currently this package is only available via Github. It is in active development and thus caution is warraranted to anyone interested in using it.
 In order to install the software use the following lines of R code

 ```r
 if(!require(devtools)){
	install.packages("devtools")
	library(devtools)
 }

install_github("apeterson91/rstapDP",dependencies = TRUE)
 ```

## Contributing

 Examples and code contributions are welcome. Feel free to start/address a feature in the issue tracker and I'll be notified shortly. 

#### Code of Conduct

Please note that `rstapDP` is released with a [Contributor Code of Conduct](https://www.contributor-covenant.org/). By contributing to this project, you agree to abide by its terms.


## How to cite this package

 A citation is in progress. Check back soon.

## Acknowledgments

This work was developed with support from NIH grant R01-HL131610 (PI: Sanchez).


Special thanks to Adam Youncan for help with the package hex.
