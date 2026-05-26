# DPTM grid-search examples
#
# This file compares the two deterministic grid-search modes for multiple
# thresholds. For a final empirical run, set grids = 100. The smaller value
# below keeps the example quick enough for a first check.

library(DPTM)

data(d1)

grids_used <- 8
grid_iter_used <- 1

# One-threshold grid search. For Th = 1, jointly and sequential are identical.
model_one <- DPTS(
  y ~ x, y ~ z,
  data = d1,
  index = c("id", "year"),
  q = d1$q,
  Th = 1,
  grid_search = TRUE,
  grids = grids_used,
  grid_search_type = "jointly"
)

print(model_one$thresholds)
print(model_one$grid_points)

# Joint two-threshold grid search.
# With grids = 100, this evaluates 100 * 100 threshold pairs.
model_joint <- DPTS(
  y ~ x, y ~ z,
  data = d1,
  index = c("id", "year"),
  q = d1$q,
  Th = 2,
  grid_search = TRUE,
  grids = grids_used,
  grid_search_type = "jointly",
  grid_search_iter = grid_iter_used
)

print(model_joint$thresholds)
print(model_joint$grid_points)

# Sequential two-threshold grid search.
# With grids = 100 and grid_search_iter = 1, this evaluates the initial
# 100 + 100 candidates plus one conditional refinement cycle.
model_seq <- DPTS(
  y ~ x, y ~ z,
  data = d1,
  index = c("id", "year"),
  q = d1$q,
  Th = 2,
  grid_search = TRUE,
  grids = grids_used,
  grid_search_type = "sequential",
  grid_search_iter = grid_iter_used
)

print(model_seq$thresholds)
print(model_seq$grid_points)

# Bootstrap test using the same grid-search interface.
# Increase bt for empirical work.
test_seq <- Threshold_Test(
  y ~ x, y ~ z,
  data = d1,
  index = c("id", "year"),
  q = d1$q,
  Th = 1,
  bt = 5,
  grid_search = TRUE,
  grids = grids_used,
  grid_search_type = "sequential",
  grid_search_iter = grid_iter_used,
  parallel = FALSE,
  seed = 2024
)

print(test_seq)
