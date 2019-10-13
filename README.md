# emba

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/bblodfon/emba.svg?branch=master)](https://travis-ci.org/bblodfon/emba)
<!-- badges: end -->

This package includes a set of functions that help the analysis and visualization of an ensemble of 
boolean models for biomarker discovery in cancer cell networks. These functions allow to easily 
import the data results of the DrugLogics software pipeline for drug combination
predictions while others (more generic ones) can be used to split a boolean model 
dataset to model groups with regards to the models predictive performance (number of true 
positive predictions or Matthews correlation coefficient score) or synergy prediction based on a given set 
of observed synergies and find the average activity difference per network 
node between all group pairs. Thus, given user-specific thresholds,
important nodes (biomarkers) can be accessed in the sense that they make the 
models predict specific synergies (synergy biomarkers) or have better 
performance in general (performance biomarkers). Lastly, if the 
boolean models have a specific equation form and differ only in their link operator, 
link operator biomarkers can also be assessed.

For an example usage of this package's functions, see the analysis performed
on multiple boolean model datasets [here](https://bblodfon.github.io/gitsbe-model-analysis/atopo/cell-lines-2500/).

## Install

```
devtools::install_github("bblodfon/emba")
```
