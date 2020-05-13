# emba

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/bblodfon/emba.svg?branch=master)](https://travis-ci.com/bblodfon/emba)
[![CRAN status](https://www.r-pkg.org/badges/version/emba)](https://cran.r-project.org/package=emba)
[![Downloads](https://cranlogs.r-pkg.org/badges/emba)](https://cran.r-project.org/package=emba)
<!-- badges: end -->

Analysis and visualization of an ensemble of boolean models for biomarker discovery in cancer cell networks. 
The package allows to easily import the data results of a software pipeline that predicts synergistic drug combinations in cancer cell lines, developed by the DrugLogics research group in NTNU. 
It has generic functions that can be used to split a boolean model 
dataset to model groups with regards to the models predictive performance (number of true 
positive predictions or Matthews correlation coefficient score) or synergy prediction based on a given set 
of observed synergies and find the average activity difference per network 
node between all group pairs. Thus, given user-specific thresholds,
important nodes (biomarkers) can be accessed in the sense that they make the 
models predict specific synergies (synergy biomarkers) or have better 
performance in general (performance biomarkers). Lastly, if the 
boolean models have a specific equation form and differ only in their link operator, 
link operator biomarkers can also be assessed.

## Install

CRAN version:
```
install.packages("emba")
```

Development version:
```
devtools::install_github("bblodfon/emba")
```

## Examples

For an example usage of this package's functions (version `0.1.1`), see [the analysis](https://bblodfon.github.io/gitsbe-model-analysis/atopo/cell-lines-2500/) performed on multiple boolean model datasets.
