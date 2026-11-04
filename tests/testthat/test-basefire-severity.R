make_severity_cache <- function(
    fuel = c(-1, 10, 20, 40, 80),
    wind = c(-1, 5, 10, 20, 40),
    last_fire_year = 0,
    last_wind_year = NULL) {

  x <- list(
    cell_id = 101L,
    region_index = 1L,
    last_fire_year = last_fire_year,
    has_wind_state = !is.null(last_wind_year),
    last_wind_year = last_wind_year,
    fuel_curve = matrix(
      as.integer(fuel),
      nrow = 1,
      dimnames = list('r', paste0('severity', 1:5))
    ),
    wind_curve = matrix(
      as.integer(wind),
      nrow = 1,
      dimnames = list('r', paste0('severity', 1:5))
    )
  )

  class(x) <- 'landisR_fire_basefire_severity_cache'
  x
}

test_that('Fuel severity increases with time since fire thresholds', {
  p <- make_severity_cache(
    fuel = c(-1, 10, 20, 40, 80)
  )

  expect_equal(.basefire_calc_severity(1, 5, p), 0L)
  expect_equal(.basefire_calc_severity(1, 10, p), 2L)
  expect_equal(.basefire_calc_severity(1, 20, p), 3L)
  expect_equal(.basefire_calc_severity(1, 40, p), 4L)
  expect_equal(.basefire_calc_severity(1, 80, p), 5L)
})

test_that('Fuel severity uses >= and later thresholds overwrite earlier ones', {
  p <- make_severity_cache(
    fuel = c(0, 0, 0, 0, 0)
  )

  expect_equal(.basefire_calc_severity(1, 0, p), 5L)
})

test_that('Minus one disables individual fuel thresholds', {
  p <- make_severity_cache(
    fuel = c(-1, -1, 30, -1, 100)
  )

  expect_equal(.basefire_calc_severity(1, 20, p), 0L)
  expect_equal(.basefire_calc_severity(1, 30, p), 3L)
  expect_equal(.basefire_calc_severity(1, 100, p), 5L)
})

test_that('Wind severity is ignored when no TimeOfLastWind exists', {
  p <- make_severity_cache(
    fuel = c(-1, -1, -1, -1, -1),
    wind = c(100, 100, 100, 100, 100),
    last_wind_year = NULL
  )

  expect_equal(.basefire_calc_severity(1, 10, p), 0L)
})

test_that('Wind severity uses <= thresholds exactly like C#', {
  p <- make_severity_cache(
    fuel = rep(-1, 5),
    wind = c(5, 10, 20, 40, 80),
    last_wind_year = 0
  )

  # TS wind = 81, no threshold passes.
  expect_equal(.basefire_calc_severity(1, 81, p), 0L)

  # TS wind = 80. S5 passes and overwrites all lower severities.
  expect_equal(.basefire_calc_severity(1, 80, p), 5L)

  # Very recent wind also passes every non-disabled threshold.
  expect_equal(.basefire_calc_severity(1, 1, p), 5L)
})

test_that('Final Base Fire severity is max of fuel and wind severity', {
  p <- make_severity_cache(
    fuel = c(-1, 10, 20, 40, 80),
    wind = c(-1, 5, 10, 15, 20),
    last_fire_year = 0,
    last_wind_year = 35
  )

  # At current 40:
  # fuel TSF=40 -> severity 4
  # wind TSW=5 -> severity 5
  out <- .basefire_calc_severity(
    1,
    40,
    p,
    diagnostics = TRUE
  )

  expect_equal(out$fuel_severity, 4L)
  expect_equal(out$wind_severity, 5L)
  expect_equal(out$severity, 5L)
})

test_that('Dynamic last-fire and last-wind vectors can override prepared state', {
  p <- make_severity_cache(
    fuel = c(-1, 10, 20, 40, 80),
    wind = c(-1, 5, 10, 20, 40),
    last_fire_year = 0,
    last_wind_year = 0
  )

  out <- .basefire_calc_severity(
    1,
    current_time = 100,
    prepared = p,
    last_fire_year = 90,
    last_wind_year = 95,
    diagnostics = TRUE
  )

  expect_equal(out$time_since_last_fire, 10)
  expect_equal(out$fuel_severity, 2L)
  expect_equal(out$time_since_last_wind, 5)
  expect_equal(out$wind_severity, 5L)
  expect_equal(out$severity, 5L)
})

test_that('Severity preparation aligns optional wind state', {
  landscape <- data.frame(
    cell_id = c(10L, 20L),
    row = c(1L, 1L),
    col = c(1L, 2L),
    active = TRUE,
    fire_region_map_code = 1L,
    max_cohort_age = c(20, 30)
  )

  topology <- fire_grid_topology(landscape)

  parameters <- .basefire_parameters(
    timestep = 10,
    species = data.frame(
      species_code = 'PICE.MAR',
      fire_tolerance = 3
    ),
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
    fuel_curves = data.frame(
      region_name = 'r',
      severity1 = -1,
      severity2 = 10,
      severity3 = 20,
      severity4 = 40,
      severity5 = 80
    ),
    wind_curves = data.frame(
      region_name = 'r',
      severity1 = -1,
      severity2 = 5,
      severity3 = 10,
      severity4 = 20,
      severity5 = 40
    ),
    damage = data.frame(
      max_age = 1,
      severity_tolerance_diff = 4
    ),
    map_names_template = 'severity.tif',
    log_file = 'events.csv',
    summary_log_file = 'summary.csv'
  )

  cache <- .basefire_prepare_parameter_cache(parameters)

  state <- fire_state(
    cell_id = c(20L, 10L),
    last_fire_year = c(-30, -20)
  )

  state$cells$last_wind_year <- c(90, 80)

  prepared <- .basefire_prepare_severity(
    landscape = landscape,
    fire_state = state,
    topology = topology,
    context = list(
      basefire_cache = cache,
      start_time = 0
    ),
    model = NULL
  )

  expect_true(prepared$has_wind_state)
  expect_equal(prepared$last_wind_year, c(80, 90))
  expect_equal(prepared$last_fire_year, c(-20, -30))
})

test_that('Severity preparation can use external TimeOfLastWind context', {
  landscape <- data.frame(
    cell_id = c(10L, 20L),
    row = c(1L, 1L),
    col = c(1L, 2L),
    active = TRUE,
    fire_region_map_code = 1L,
    max_cohort_age = c(20, 30)
  )

  topology <- fire_grid_topology(landscape)

  parameters <- .basefire_parameters(
    timestep = 10,
    species = data.frame(
      species_code = 'PICE.MAR',
      fire_tolerance = 3
    ),
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

  cache <- .basefire_prepare_parameter_cache(parameters)

  state <- fire_state(
    cell_id = c(10L, 20L),
    last_fire_year = c(-20, -30)
  )

  prepared <- .basefire_prepare_severity(
    landscape = landscape,
    fire_state = state,
    topology = topology,
    context = list(
      basefire_cache = cache,
      start_time = 0,
      time_of_last_wind = c('10' = 70, '20' = 90)
    ),
    model = NULL
  )

  expect_true(prepared$has_wind_state)
  expect_equal(prepared$last_wind_year, c(70, 90))
})

test_that('Base Fire severity is a real internal severity component', {
  expect_true(is.function(.basefire_severity))
  expect_equal(
    attr(.basefire_severity, 'fire_component_class'),
    'severity'
  )
  expect_true(is.function(
    attr(.basefire_severity, 'fire_prepare')
  ))
})
