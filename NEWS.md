# emba 0.1.5

- Finally added tests to the package! **Coverage is now 76%**.
- Used the [pkgdown](https://github.com/r-lib/pkgdown/) package to create static html documentation for emba. [Check it here](https://bblodfon.github.io/emba/index.html)!
- Change MCC calculation to return 0 when undefined/`NaN` MCC scores were produced (which is the correct limiting value - see [Chicco at al. (2020)](https://doi.org/10.1186/s12864-019-6413-7)). Thus, the previous versions handling of `NaN` MCC scores, is now deprecated.
- Add the `penalty` parameter to account for the difference in model group size when calculating the average activity or link operator data differences. This minimizes the bias in the returned biomarkers.
    - For the implementation check the function `emba::get_vector_diff` and the corresponding [StackOverflow question](https://math.stackexchange.com/questions/3547139/formula-for-weighted-average-difference).
    - To get the same results as with previous versions of this library, use `penalty=0` in the general `emba::biomarker_*` functions (though the results will probably be very biased and that's why the default value for the `penalty` is now **0.1**).
- The following functions do not take the input parameter `models` anymore:
  - `emba::get_avg_link_operator_diff_mat_based_on_tp_predictions`
  - `emba::get_avg_activity_diff_mat_based_on_tp_predictions`
  - `emba::get_avg_activity_diff_based_on_tp_predictions`
- Refactor several of the functions that load the results from the DrugLogics pipeline:
  - If a model has less or more than 1 stable state, it's discarded and a message is printed.
  - Return value is now a `data.frame` object instead of a `matrix`.
  - The models names do not have the annoying `.gitsbe` extension anymore.
  - These changes affect the following functions: `emba::get_link_operators_from_models_dir`, `emba::get_stable_state_from_models_dir` and `emba::get_model_names`.
- Plot (`emba::plot_*`) functions refactoring
- Various small bug fixes :)

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
