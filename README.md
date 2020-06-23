# emba

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/bblodfon/emba.svg?branch=master)](https://travis-ci.com/bblodfon/emba)
[![codecov](https://codecov.io/gh/bblodfon/emba/branch/master/graph/badge.svg)](https://codecov.io/gh/bblodfon/emba)
[![CRAN status](https://www.r-pkg.org/badges/version/emba)](https://cran.r-project.org/package=emba)
[![Downloads](https://cranlogs.r-pkg.org/badges/emba)](https://cran.r-project.org/package=emba)
<!-- badges: end -->

Analysis and visualization of an ensemble of boolean models for biomarker discovery in cancer cell networks. 

The package allows to easily load the simulation data results of the [DrugLogics](https://github.com/druglogics) software pipeline that is used to 
predict synergistic drug combinations in cancer cell lines. 
It has generic functions that can be used to split a boolean model 
dataset to model groups with regards to the models predictive performance (number of *true 
positive* predictions/*Matthews correlation coefficient* score) or synergy prediction based on a given set 
of observed synergies and find the average activity difference per network 
node between all model group pairs. Thus, given user-specific thresholds,
important nodes (*biomarkers*) can be accessed in the sense that they make the 
models predict specific synergies (*synergy biomarkers*) or have better 
performance in general (*performance biomarkers*). 

Lastly, if the boolean models have a [specific equation form](https://druglogics.github.io/druglogics-doc/gitsbe-description.html#default-equation) and differ only in their link operator, *link operator* biomarkers can also be found.

## Install

CRAN version:
```
install.packages("emba")
```

Development version:
```
devtools::install_github("bblodfon/emba")
```

## Usage

Check the [Get Started guide](https://bblodfon.github.io/emba/articles/emba.html).

For an earlier example usage of this package (version `0.1.1`), see this [analysis](https://bblodfon.github.io/gitsbe-model-analysis/atopo/cell-lines-2500/) performed on multiple boolean model datasets.

## Code of Conduct

Please note that the emba project is released with a [Contributor Code of Conduct](https://bblodfon.github.io/emba/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
