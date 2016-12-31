
# grpregOverlap
[![Build Status](https://travis-ci.org/YaohuiZeng/grpregOverlap.svg?branch=master)](https://travis-ci.org/YaohuiZeng/grpregOverlap)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/grpregOverlap)](https://CRAN.R-project.org/package=grpregOverlap)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/grpregOverlap)](http://www.r-pkg.org/pkg/grpregOverlap)

`grpregOverlap` fits the regularization paths of linear, logistic, Poisson, 
or Cox models with overlapping grouped covariates based on the latent group lasso 
approach (Jacob et al., 2009). Latent group MCP/SCAD as well as bi-level 
selection methods, namely the group exponential lasso(Breheny, 2015) and the 
composite MCP (Huang et al., 2012) are also available. This package serves as 
an extension of R package `grpreg` (by Dr. Patrick Breheny <patrick-breheny@uiowa.edu>) 
for grouped variable selection involving overlaps between groups.

## News:
* this package now works for survival analysis (Cox model) by specifying "family = cox". 
* this package on GitHub has been updated to Version 2.2-0. See details in NEWS.

## Installation:
* the stable version: `install.packages("grpregOverlap")`
* the latest version: `devtools::install_github("YaohuiZeng/grpregOverlap")`

## Report bugs：
* open an [issue](https://github.com/YaohuiZeng/grpregOverlap/issues) or send an email to Yaohui Zeng at <yaohui-zeng@uiowa.edu>
