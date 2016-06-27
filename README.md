[![Build Status](https://travis-ci.org/YaohuiZeng/grpregOverlap.svg?branch=master)](https://travis-ci.org/YaohuiZeng/grpregOverlap)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/grpregOverlap)](http://cran.r-project.org/package=grpregOverlap)

# grpregOverlap

`grpregOverlap` fits the regularization paths of linear, logistic, Poisson, 
or Cox models with overlapping grouped covariates based on the latent group lasso 
approach (Jacob et al., 2009). Latent group MCP/SCAD as well as bi-level 
selection methods, namely the group exponential lasso(Breheny, 2015) and the 
composite MCP (Huang et al., 2012) are also available. This package serves as 
an extension of R package `grpreg` (by Dr. Patrick Breheny <patrick-breheny@uiowa.edu>) 
for grouped variable selection involving overlaps between groups.

To install:
* the latest version (requires devtools): install_github("YaohuiZeng/grpregOverlap")

To report bugs：
* send email to Yaohui Zeng at <yaohui-zeng@uiowa.edu>
