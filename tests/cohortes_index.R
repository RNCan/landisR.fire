kernel <- getFromNamespace("basefire_timestep_kernel_260", "landisR.fire")
grid <- expand.grid(row = 1:3, col = 1:3)
grid$cell_id <- seq_len(nrow(grid)) * 10L
topology <- landisR.fire::fire_grid_topology(grid)
n <- nrow(grid)

run <- function(cohort_sites, use_offsets) {
  set.seed(4357)
  cohort_env <- new.env(parent = emptyenv())
  cohort_env$alive <- rep(TRUE, length(cohort_sites))
  index <- cohort_sites
  if (!use_offsets) {
    groups <- split(seq_along(cohort_sites), cohort_sites)
    index <- unname(groups[as.character(seq_len(n))])
  }

  fire <- kernel(
    current_time = 10L, timestep = 10,
    ignition_probability = rep(0.1, n),
    ignition_fire_spread_age = rep(5, n),
    cell_id = grid$cell_id, region_index = rep(1L, n), n_regions = 1L,
    neighbors_full = topology$neighbors,
    second_neighbors_full = topology$second_neighbors,
    full_to_active = seq_len(n), full_active = rep(TRUE, n),
    fire_spread_age = rep(5, n), fuel_curve = matrix(0, 1L, 5L),
    wind_curve = matrix(-1, 1L, 5L), has_wind = FALSE, last_wind_year = numeric(),
    last_fire_year_input = rep(-50, n), event_occupied_input = rep(FALSE, n),
    last_severity_input = rep(NA_real_, n), rows_by_site = index,
    cohort_env = cohort_env, cohort_age = rep(50, length(cohort_sites)),
    cohort_species_index = rep(1L, length(cohort_sites)),
    species_fire_tolerance = 1L, species_longevity = 100, damage_max_age = 1,
    damage_severity_tolerance_diff = 0L, size_one = function(region) n
  )
  list(fire = fire, alive = cohort_env$alive, rng = .Random.seed)
}

ordered <- c(1L, 1L, 3L, 3L, 4L, 7L, 9L)
interleaved <- c(9L, 3L, 4L, 1L, 3L, 7L, 1L, 0L, 9L)
for (sites in list(ordered, interleaved, integer(), c(0L, 0L))) {
  reference <- run(sites, use_offsets = FALSE)
  candidate <- run(sites, use_offsets = TRUE)
  stopifnot(identical(candidate, reference), candidate$fire$number_fires > 0L,
            identical(candidate$alive, sites == 0L))
}
