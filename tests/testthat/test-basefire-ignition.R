make_uniform_sequence <- function(values) {
  i <- 0L

  function(n = 1L) {
    if (!identical(as.integer(n), 1L)) {
      stop('test uniform supports n=1 only')
    }

    i <<- i + 1L

    if (i > length(values)) {
      stop('uniform sequence exhausted')
    }

    values[i]
  }
}

make_basefire_ignition_cache <- function(
    ignition_probability = 0.1,
    fire_spread_age = 20,
    last_fire_year = 0) {

  x <- list(
    cell_id = 101L,
    region_index = 1L,
    ignition_probability = ignition_probability,
    fire_spread_age = fire_spread_age,
    last_fire_year = last_fire_year
  )

  class(x) <- 'landisR_fire_basefire_ignition_cache'
  x
}

test_that('Base Fire initiation/spread probability matches reference formula', {
  expect_equal(
    .basefire_compute_init_spread_prob(20, 20),
    1 - exp(-1)
  )

  expect_equal(
    .basefire_compute_init_spread_prob(0, 20),
    0
  )

  expect_equal(
    .basefire_compute_init_spread_prob(10, 0),
    1
  )

  expect_true(
    is.nan(.basefire_compute_init_spread_prob(0, 0))
  )
})

test_that('Base Fire initializes TimeOfLastFire from maximum cohort age', {
  expect_equal(
    .basefire_initial_last_fire_year(
      max_cohort_age = c(0, 25, 80),
      start_time = 0
    ),
    c(0, -25, -80)
  )

  expect_equal(
    .basefire_initial_last_fire_year(
      max_cohort_age = c(10, 30),
      start_time = 100
    ),
    c(90, 70)
  )
})

test_that('Ignition gate failure consumes exactly one RNG draw', {
  prepared <- make_basefire_ignition_cache(
    ignition_probability = 0.1,
    fire_spread_age = 20,
    last_fire_year = 0
  )

  severity_calls <- 0L
  severity_fun <- function(cell_index, current_time) {
    severity_calls <<- severity_calls + 1L
    3
  }

  out <- .basefire_try_initiate(
    cell_index = 1,
    current_time = 20,
    timestep = 1,
    prepared = prepared,
    severity_fun = severity_fun,
    uniform = make_uniform_sequence(c(0.5, 0.1))
  )

  expect_false(out$initiated)
  expect_equal(out$failed_gate, 'ignition')
  expect_equal(out$draws, 1L)
  expect_equal(severity_calls, 0L)
})

test_that('Initiation/spread gate failure consumes two RNG draws', {
  prepared <- make_basefire_ignition_cache(
    ignition_probability = 0.5,
    fire_spread_age = 20,
    last_fire_year = 0
  )

  severity_calls <- 0L
  severity_fun <- function(cell_index, current_time) {
    severity_calls <<- severity_calls + 1L
    3
  }

  out <- .basefire_try_initiate(
    cell_index = 1,
    current_time = 20,
    timestep = 1,
    prepared = prepared,
    severity_fun = severity_fun,
    uniform = make_uniform_sequence(c(0.1, 0.9))
  )

  expect_false(out$initiated)
  expect_equal(out$failed_gate, 'init_spread')
  expect_equal(out$draws, 2L)
  expect_equal(severity_calls, 0L)
})

test_that('Severity is evaluated only after both probability gates', {
  prepared <- make_basefire_ignition_cache(
    ignition_probability = 0.5,
    fire_spread_age = 20,
    last_fire_year = 0
  )

  severity_calls <- 0L
  severity_fun <- function(cell_index, current_time) {
    severity_calls <<- severity_calls + 1L
    0
  }

  out <- .basefire_try_initiate(
    cell_index = 1,
    current_time = 20,
    timestep = 1,
    prepared = prepared,
    severity_fun = severity_fun,
    uniform = make_uniform_sequence(c(0.1, 0.1))
  )

  expect_false(out$initiated)
  expect_equal(out$failed_gate, 'severity')
  expect_equal(out$draws, 2L)
  expect_equal(severity_calls, 1L)
  expect_equal(out$severity, 0)
})

test_that('Base Fire initiation succeeds when all three gates pass', {
  prepared <- make_basefire_ignition_cache(
    ignition_probability = 0.5,
    fire_spread_age = 20,
    last_fire_year = 0
  )

  out <- .basefire_try_initiate(
    cell_index = 1,
    current_time = 20,
    timestep = 1,
    prepared = prepared,
    severity_fun = function(cell_index, current_time) 2,
    uniform = make_uniform_sequence(c(0.1, 0.1))
  )

  expect_true(out$initiated)
  expect_true(is.na(out$failed_gate))
  expect_equal(out$draws, 2L)
  expect_equal(out$severity, 2)
})

test_that('IgnitionProbability multiplied by timestep is not clamped', {
  prepared <- make_basefire_ignition_cache(
    ignition_probability = 0.2,
    fire_spread_age = 20,
    last_fire_year = 0
  )

  out <- .basefire_try_initiate(
    cell_index = 1,
    current_time = 20,
    timestep = 10,
    prepared = prepared,
    severity_fun = function(cell_index, current_time) 2,
    uniform = make_uniform_sequence(c(0.999, 0.1))
  )

  expect_equal(out$ignition_probability, 2)
  expect_true(out$initiated)
})

test_that('Base Fire uses <= boundary for both probability gates', {
  prepared <- make_basefire_ignition_cache(
    ignition_probability = 0.25,
    fire_spread_age = 20,
    last_fire_year = 0
  )

  p2 <- .basefire_compute_init_spread_prob(20, 20)

  out <- .basefire_try_initiate(
    cell_index = 1,
    current_time = 20,
    timestep = 1,
    prepared = prepared,
    severity_fun = function(cell_index, current_time) 1,
    uniform = make_uniform_sequence(c(0.25, p2))
  )

  expect_true(out$initiated)
})

test_that('FireSpreadAge zero at TSF zero cannot pass NaN gate', {
  prepared <- make_basefire_ignition_cache(
    ignition_probability = 1,
    fire_spread_age = 0,
    last_fire_year = 20
  )

  out <- .basefire_try_initiate(
    cell_index = 1,
    current_time = 20,
    timestep = 1,
    prepared = prepared,
    severity_fun = function(cell_index, current_time) 5,
    uniform = make_uniform_sequence(c(0.1, 0))
  )

  expect_true(is.nan(out$init_spread_probability))
  expect_false(out$initiated)
  expect_equal(out$failed_gate, 'init_spread')
})

test_that('Base Fire ignition preparation aligns state and region map', {
  landscape <- data.frame(
    cell_id = c(10L, 20L, 30L),
    row = c(1L, 1L, 1L),
    col = c(1L, 2L, 3L),
    active = TRUE,
    fire_region_map_code = c(2L, 5L, 2L),
    max_cohort_age = c(10, 20, 30)
  )

  regions <- .basefire_normalize_regions(data.frame(
    name = c('r2', 'r5'),
    map_code = c(2, 5),
    mean_size = c(100, 200),
    min_size = c(1, 1),
    max_size = c(1000, 2000),
    ignition_probability = c(0.01, 0.02),
    fire_spread_age = c(20, 40)
  ))

  p <- .basefire_parameters(
    timestep = 10,
    species = data.frame(
      species_code = 'PICE.MAR',
      fire_tolerance = 3
    ),
    regions = regions[, c(
      'name',
      'map_code',
      'mean_size',
      'min_size',
      'max_size',
      'ignition_probability',
      'fire_spread_age'
    )],
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

  base_cache <- .basefire_prepare_parameter_cache(p)

  state <- fire_state(
    cell_id = c(30L, 10L, 20L),
    last_fire_year = c(NA, -99, NA)
  )

  topology <- fire_grid_topology(landscape)

  prepared <- .basefire_prepare_ignition(
    landscape = landscape,
    fire_state = state,
    topology = topology,
    context = list(
      basefire_cache = base_cache,
      start_time = 0
    ),
    model = NULL
  )

  expect_equal(prepared$cell_id, c(10L, 20L, 30L))
  expect_equal(prepared$region_index, c(1L, 2L, 1L))
  expect_equal(prepared$ignition_probability, c(0.01, 0.02, 0.01))
  expect_equal(prepared$fire_spread_age, c(20L, 40L, 20L))

  # Cell 10 retains explicit fire state; missing values use StartTime - max age.
  expect_equal(prepared$last_fire_year, c(-99, -20, -30))
})

test_that('Base Fire ignition is now a real internal component', {
  expect_true(is.function(.basefire_ignition))
  expect_equal(
    attr(.basefire_ignition, 'fire_component_class'),
    'ignition'
  )
  expect_true(is.function(
    attr(.basefire_ignition, 'fire_prepare')
  ))
})
