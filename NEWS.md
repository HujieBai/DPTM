# DPTM 4.0.1

* Added `grid_search_iter` to `DPTS()` and `Threshold_Test()` for coordinate-refinement cycles after deterministic grid search.
* For multiple thresholds, one refinement cycle re-estimates each threshold conditional on the other current thresholds and sorts the final thresholds.
* Printed grid-search results now report the number of refinement iterations.

# DPTM 4.0.0

* Added `grid_search`, `grids`, and `grid_search_type` to `DPTS()` and `Threshold_Test()` so unknown thresholds can be searched by MCMC, joint grid search, or sequential grid search.
* Shared threshold-design construction across fixed-threshold fitting, MCMC search, and grid search to reduce duplicated model-building code.
* Improved bootstrap retry handling in `Threshold_Test()` and made the parallel backend shutdown explicit.
* Clarified README and help pages for MCMC versus grid-search threshold estimation.

# DPTM 3.0.2

# DPTM 3.0.0

* Rewrote the entire package in OOP (R6 class).
* Detailed instructions and examples have been added to make it more user-friendly.
* The formula interface of R language is exploited, which is more standardized and easy to use.

# DPTM 2.0.0
