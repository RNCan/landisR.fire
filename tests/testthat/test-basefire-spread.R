make_uniform_counter <- function(value = 0) {
  e <- new.env(parent = emptyenv())
  e$n <- 0L

  fun <- function(n = 1L) {
    e$n <- e$n + 1L
    rep(value, n)
  }

  attr(fun, 'counter') <- e
  fun
}

make_sequence_draw <- function(values) {
  i <- 0L

  function(n = 1L) {
    if (n != 1L) stop('test draw supports n=1 only')

    i <<- i + 1L

    if (i > length(values)) {
      stop('draw sequence exhausted')
    }

    values[i]
  }
}

make_spread_fixture <- function(inactive_cell = integer()) {
  grid <- expand.grid(
    row = 1:5,
    col = 1:5
  )

  grid <- grid[order(grid$row, grid$col), , drop = FALSE]
  grid$cell_id <- seq_len(nrow(grid))
  grid$active <- !(grid$cell_id %in% inactive_cell)
  grid$fire_region_map_code <- 1L
  grid$max_cohort_age <- 20

  topology <- fire_grid_topology(grid)

  parameters <- .basefire_parameters(
    timestep = 10,
    species = data.frame(
      species_code = 'PICE.MAR',
      fire_tolerance = 3
    ),
    regions = data.frame(
      name = 'region-1',
      map_code = 1,
      mean_size = 100,
      min_size = 1,
      max_size = 1000,
      ignition_probability = 0.01,
      fire_spread_age = 20
    ),
    initial_fire_regions = 'regions.tif',
    fuel_curves = NULL,
    wind_curves = NULL,
    damage = data.frame(
      max_age = 1,
      severity_tolerance_diff = 4
    ),
    map_names_template = 'severity.tif',
    log_file = 'events.csv',
    summary_log_file = 'summary.csv'
  )

  base_cache <- .basefire_prepare_parameter_cache(parameters)

  spread <- .basefire_prepare_spread(
    landscape = grid,
    fire_state = NULL,
    topology = topology,
    context = list(
      basefire_cache = base_cache,
      basefire_full_landscape = grid
    ),
    model = NULL
  )

  burn <- .basefire_prepare_burn_probability(
    landscape = grid,
    fire_state = NULL,
    topology = topology,
    context = list(
      basefire_cache = base_cache,
      start_time = 0
    ),
    model = NULL
  )

  size <- .basefire_prepare_size(
    landscape = grid,
    fire_state = NULL,
    topology = topology,
    context = list(
      basefire_cache = base_cache,
      cell_area = 1
    ),
    model = NULL
  )

  list(
    grid = grid,
    topology = topology,
    base_cache = base_cache,
    spread = spread,
    burn = burn,
    size = size
  )
}

test_that('Base Fire wind probabilities reproduce GetNeighbors formula', {
  expect_equal(
    .basefire_wind_probabilities(0),
    rep(0.5, 8)
  )

  expect_equal(
    .basefire_wind_probabilities(1),
    c(0.75, 0.75, 0.375, 0, 0, 0, 0.375, 0.75)
  )
})

test_that('Base Fire site topology distinguishes null and inactive sites', {
  grid <- expand.grid(row = 1:3, col = 1:3)
  grid <- grid[order(grid$row, grid$col), , drop = FALSE]
  grid$cell_id <- seq_len(9)
  grid$active <- TRUE

  east_id <- grid$cell_id[grid$row == 2 & grid$col == 3]
  grid$active[grid$cell_id == east_id] <- FALSE

  active_topology <- fire_grid_topology(grid)

  x <- .basefire_build_site_topology(
    landscape = grid,
    active_topology = active_topology,
    full_landscape = grid
  )

  center_id <- grid$cell_id[grid$row == 2 & grid$col == 2]
  center_active <- match(center_id, active_topology$cell_id)

  east_full <- x$neighbors_full[center_active, 'E']

  expect_gt(east_full, 0L)
  expect_false(x$full_active[east_full])
  expect_equal(x$full_to_active[east_full], 0L)
})

test_that('GetNeighbors consumes draws for inactive sites but not null sites', {
  grid <- expand.grid(row = 1:3, col = 1:3)
  grid <- grid[order(grid$row, grid$col), , drop = FALSE]
  grid$cell_id <- seq_len(9)
  grid$active <- TRUE

  east_id <- grid$cell_id[grid$row == 2 & grid$col == 3]
  grid$active[grid$cell_id == east_id] <- FALSE
  grid$fire_region_map_code <- 1L

  topology <- fire_grid_topology(grid)

  parameters <- .basefire_parameters(
    timestep = 10,
    species = data.frame(species_code = 'PICE.MAR', fire_tolerance = 3),
    regions = data.frame(
      name = 'r',
      map_code = 1,
      mean_size = 100,
      min_size = 1,
      max_size = 1000,
      ignition_probability = 0.01,
      fire_spread_age = 20
    ),
    initial_fire_regions = 'regions.tif',
    damage = data.frame(max_age = 1, severity_tolerance_diff = 4),
    map_names_template = 'severity.tif',
    log_file = 'events.csv',
    summary_log_file = 'summary.csv'
  )

  cache <- .basefire_prepare_parameter_cache(parameters)

  prepared <- .basefire_prepare_spread(
    landscape = grid,
    fire_state = NULL,
    topology = topology,
    context = list(
      basefire_cache = cache,
      basefire_full_landscape = grid
    ),
    model = NULL
  )

  center_id <- grid$cell_id[grid$row == 2 & grid$col == 2]
  center <- match(center_id, topology$cell_id)

  nd <- make_uniform_counter(0)
  gu <- make_uniform_counter(0.5)

  neighbors <- .basefire_get_neighbors(
    cell_index = center,
    wind_direction = 0,
    wind_speed = 0,
    prepared = prepared,
    next_double = nd,
    generate_uniform = gu
  )

  # Every first-ring geometric site exists, including the inactive east cell.
  expect_equal(attr(neighbors, 'next_double_draws'), 8L)
  expect_equal(nrow(neighbors), 8L)
  expect_true(any(!neighbors$is_active))
  expect_equal(attr(nd, 'counter')$n, 8L)

  # The ninth N neighbor is outside a 3x3 grid, so no GenerateUniform is used.
  expect_equal(attr(neighbors, 'uniform_draws'), 0L)
  expect_equal(attr(gu, 'counter')$n, 0L)
})

test_that('GetNeighbors adds ninth site two cells downwind', {
  fx <- make_spread_fixture()

  source_id <- fx$grid$cell_id[
    fx$grid$row == 3 & fx$grid$col == 3
  ]
  source <- match(source_id, fx$topology$cell_id)

  neighbors <- .basefire_get_neighbors(
    cell_index = source,
    wind_direction = 0,
    wind_speed = 1,
    prepared = fx$spread,
    next_double = function(n = 1L) rep(0, n),
    generate_uniform = function(n = 1L) rep(0, n)
  )

  ninth <- neighbors[neighbors$ring == 2L, , drop = FALSE]

  expect_equal(nrow(ninth), 1L)

  full_row <- fx$grid$row[
    match(ninth$cell_id, fx$grid$cell_id)
  ]
  full_col <- fx$grid$col[
    match(ninth$cell_id, fx$grid$cell_id)
  ]

  expect_equal(full_row, 1L)
  expect_equal(full_col, 3L)
  expect_equal(ninth$direction, 'N')
})

test_that('Wind direction above 7 is clamped to northwest', {
  fx <- make_spread_fixture()

  source_id <- fx$grid$cell_id[
    fx$grid$row == 3 & fx$grid$col == 3
  ]
  source <- match(source_id, fx$topology$cell_id)

  neighbors <- .basefire_get_neighbors(
    cell_index = source,
    wind_direction = 8,
    wind_speed = 1,
    prepared = fx$spread,
    next_double = function(n = 1L) rep(1, n),
    generate_uniform = function(n = 1L) rep(0, n)
  )

  ninth <- neighbors[neighbors$ring == 2L, , drop = FALSE]

  expect_equal(attr(neighbors, 'wind_direction'), 7L)
  expect_equal(nrow(ninth), 1L)
  expect_equal(ninth$direction, 'NW')
})

test_that('Burn probability is the Base Fire TSF gate', {
  fx <- make_spread_fixture()

  source_id <- fx$grid$cell_id[
    fx$grid$row == 3 & fx$grid$col == 3
  ]
  source <- match(source_id, fx$topology$cell_id)

  p <- .basefire_compute_burn_probability(
    cell_index = source,
    current_time = 20,
    prepared = fx$burn,
    last_fire_year = rep(0, length(fx$burn$cell_id))
  )

  expect_equal(p, 1 - exp(-1))
})

test_that('Spread event follows FIFO queue and region target', {
  fx <- make_spread_fixture()

  source_id <- fx$grid$cell_id[
    fx$grid$row == 3 & fx$grid$col == 3
  ]
  source <- match(source_id, fx$topology$cell_id)

  fixed_size <- function(
      region_index,
      prepared,
      exponential_draw,
      diagnostics = FALSE) {
    3L
  }

  always_burn <- function(
      cell_index,
      current_time,
      prepared,
      last_fire_year = NULL,
      ...) {
    1
  }

  identity_shuffle <- function(x) x

  out <- .basefire_run_spread_event(
    initiation_cell_index = source,
    current_time = 20,
    spread_prepared = fx$spread,
    size_prepared = fx$size,
    burn_prepared = fx$burn,
    severity_fun = function(cell_index, current_time) 1,
    damage_fun = function(cell_index, severity, current_time) 1L,
    size_fun = fixed_size,
    burn_probability_fun = always_burn,
    next_double = function(n = 1L) rep(0, n),
    generate_uniform = function(n = 1L) rep(0.5, n),
    shuffle_fun = identity_shuffle,
    wind_direction = 0,
    wind_speed = 0
  )

  expect_equal(out$total_sites_in_event, 3L)
  expect_equal(out$num_site_checked, 3L)
  expect_equal(out$sites_in_event[1], 3L)
  expect_equal(out$size_targets[1], 3L)
  expect_equal(nrow(out$cells), 3L)
  expect_true(all(out$cells$burned))
  expect_equal(out$total_sites_damaged, 3L)
  expect_equal(out$severity, 1)
})

test_that('Severity zero marks event membership but does not spread', {
  fx <- make_spread_fixture()

  source_id <- fx$grid$cell_id[
    fx$grid$row == 3 & fx$grid$col == 3
  ]
  source <- match(source_id, fx$topology$cell_id)

  out <- .basefire_run_spread_event(
    initiation_cell_index = source,
    current_time = 20,
    spread_prepared = fx$spread,
    size_prepared = fx$size,
    burn_prepared = fx$burn,
    severity_fun = function(cell_index, current_time) 0,
    size_fun = function(
      region_index,
      prepared,
      exponential_draw,
      diagnostics = FALSE) 10L,
    wind_direction = 0,
    wind_speed = 0
  )

  expect_equal(out$num_site_checked, 1L)
  expect_equal(out$total_sites_in_event, 0L)
  expect_equal(nrow(out$cells), 1L)
  expect_false(out$cells$burned)
  expect_true(out$event_occupied[source])
})

test_that('Inactive selected neighbors are shuffled but never queued', {
  grid <- expand.grid(row = 1:5, col = 1:5)
  grid <- grid[order(grid$row, grid$col), , drop = FALSE]
  grid$cell_id <- seq_len(nrow(grid))
  grid$active <- TRUE

  source_id <- grid$cell_id[grid$row == 3 & grid$col == 3]
  east_id <- grid$cell_id[grid$row == 3 & grid$col == 4]
  grid$active[grid$cell_id == east_id] <- FALSE
  grid$fire_region_map_code <- 1L
  grid$max_cohort_age <- 20

  topology <- fire_grid_topology(grid)

  parameters <- .basefire_parameters(
    timestep = 10,
    species = data.frame(species_code = 'PICE.MAR', fire_tolerance = 3),
    regions = data.frame(
      name = 'r',
      map_code = 1,
      mean_size = 100,
      min_size = 1,
      max_size = 1000,
      ignition_probability = 0.01,
      fire_spread_age = 20
    ),
    initial_fire_regions = 'regions.tif',
    damage = data.frame(max_age = 1, severity_tolerance_diff = 4),
    map_names_template = 'severity.tif',
    log_file = 'events.csv',
    summary_log_file = 'summary.csv'
  )

  cache <- .basefire_prepare_parameter_cache(parameters)

  spread <- .basefire_prepare_spread(
    landscape = grid,
    fire_state = NULL,
    topology = topology,
    context = list(
      basefire_cache = cache,
      basefire_full_landscape = grid
    ),
    model = NULL
  )

  burn <- .basefire_prepare_burn_probability(
    landscape = grid,
    fire_state = NULL,
    topology = topology,
    context = list(
      basefire_cache = cache,
      start_time = 0
    ),
    model = NULL
  )

  size <- .basefire_prepare_size(
    landscape = grid,
    fire_state = NULL,
    topology = topology,
    context = list(basefire_cache = cache, cell_area = 1),
    model = NULL
  )

  source <- match(source_id, topology$cell_id)

  seen_inactive <- FALSE

  shuffle_probe <- function(x) {
    seen_inactive <<- any(!x$is_active)
    x
  }

  out <- .basefire_run_spread_event(
    initiation_cell_index = source,
    current_time = 20,
    spread_prepared = spread,
    size_prepared = size,
    burn_prepared = burn,
    severity_fun = function(cell_index, current_time) 1,
    damage_fun = function(cell_index, severity, current_time) 1L,
    size_fun = function(
      region_index,
      prepared,
      exponential_draw,
      diagnostics = FALSE) 2L,
    burn_probability_fun = function(
      cell_index,
      current_time,
      prepared,
      last_fire_year = NULL,
      ...) 1,
    next_double = function(n = 1L) rep(0, n),
    generate_uniform = function(n = 1L) rep(0.5, n),
    shuffle_fun = shuffle_probe,
    wind_direction = 0,
    wind_speed = 0
  )

  expect_true(seen_inactive)
  expect_false(east_id %in% out$cells$cell_id)
})

test_that('Base Fire spread and burn probability are real components', {
  expect_true(is.function(.basefire_spread))
  expect_true(is.function(.basefire_burn_probability))

  expect_equal(
    attr(.basefire_spread, 'fire_component_class'),
    'spread'
  )

  expect_equal(
    attr(.basefire_burn_probability, 'fire_component_class'),
    'burn_probability'
  )

  expect_true(is.function(attr(.basefire_spread, 'fire_prepare')))
  expect_true(is.function(
    attr(.basefire_burn_probability, 'fire_prepare')
  ))
})

test_that('Occupied initiation site is not rejected by strict Base Fire spread', {
  fx <- make_spread_fixture()

  source_id <- fx$grid$cell_id[
    fx$grid$row == 3 & fx$grid$col == 3
  ]
  source <- match(source_id, fx$topology$cell_id)

  occupied <- rep(FALSE, length(fx$topology$cell_id))
  occupied[source] <- TRUE

  out <- .basefire_run_spread_event(
    initiation_cell_index = source,
    current_time = 20,
    spread_prepared = fx$spread,
    size_prepared = fx$size,
    burn_prepared = fx$burn,
    severity_fun = function(cell_index, current_time) 1,
    damage_fun = function(cell_index, severity, current_time) 0L,
    size_fun = function(
      region_index,
      prepared,
      exponential_draw,
      diagnostics = FALSE) 1L,
    event_occupied = occupied,
    next_double = function(n = 1L) rep(1, n),
    generate_uniform = function(n = 1L) rep(0.5, n),
    shuffle_fun = function(x) x,
    wind_direction = 0,
    wind_speed = 0
  )

  expect_equal(out$num_site_checked, 1L)
  expect_equal(out$total_sites_in_event, 1L)
  expect_true(out$event_occupied[source])
  expect_true(out$event_member[source])
  expect_true(out$cells$burned[out$cells$cell_index == source])
})
