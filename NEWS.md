# emba 0.1.4

- `get_synergy_scores` now supports reading both *ensemble-wise* and *model-wise* synergies files
- add `calculate.subsets.stats` option to the general analysis functions (`biomarker_*`) that decides if the powerset of the observed synergies and the number of models predicting each subset is going to be calculated. 
The default value is set to `FALSE` to save computation time :)

# emba 0.1.3

- add function `get_synergy_scores`
- fixed test that used a randomly generated matrix

# emba 0.1.2

- add function `get_avg_link_operator_diff_based_on_synergy_set_cmp`
- add function `get_avg_link_operator_diff_based_on_specific_synergy_prediction`
- add function `filter_network` - to use for visualizing induced subgraphs
- update dependencies (set `usefun` min version to 0.4.3)

# emba 0.1.1

- Optimized `count_models_that_predict_synergies` function and added tests for it. For a benchmark see
relative [Stack Overflow thread](https://stackoverflow.com/questions/58380043/optimize-r-code-for-row-operations-on-ternary-data-frame)

# emba 0.1.0

- Added a `NEWS.md` file to track changes to the package
- Transferred functions from separate R scripts to the package
- Finished code refactoring and splitting to different modules
- Finished writing documentation for all functions
