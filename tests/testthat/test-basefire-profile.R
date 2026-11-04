make_basefire_profile_fixture <- function() {
  landscape <- expand.grid(
    row = 1:3,
    col = 1:3
  )

  landscape <- landscape[
    order(landscape$row, landscape$col),
    ,
    drop = FALSE
  ]

  landscape$cell_id <- seq_len(nrow(landscape))
  landscape$active <- TRUE
  landscape$fire_region_map_code <- 1L
  landscape$max_cohort_age <- 50

  parameters <- .basefire_parameters(
    timestep = 10,
    species = data.frame(
      species_code = 'PICE.MAR',
      fire_tolerance = 3
    ),
    regions = data.frame(
      name = 'r',
      map_code = 1,
      mean_size = 3,
      min_size = 3,
      max_size = 3,
      ignition_probability = 0.01,
      fire_spread_age = 20
    ),
    initial_fire_regions = 'regions.tif',
    fuel_curves = data.frame(
      region_name = 'r',
      severity1 = 0,
      severity2 = 20,
      severity3 = 100,
      severity4 = 150,
      severity5 = 200
    ),
    wind_curves = NULL,
    damage = data.frame(
      max_age = c(0.25, 1),
      severity_tolerance_diff = c(-2, 4)
    ),
    map_names_template = 'severity-{timestep}.tif',
    log_file = 'events.csv',
    summary_log_file = 'summary.csv'
  )

  profile <- .basefire_profile(parameters)

  prepared <- .basefire_prepare_profile(
    profile = profile,
    landscape = landscape,
    fire_state = NULL,
    context = list(
      start_time = 0,
      cell_area = 1,
      species_longevity = c(
        'PICE.MAR' = 100
      )
    )
  )

  cohorts <- data.frame(
    cell_id = rep(landscape$cell_id, each = 2),
    species_code = 'PICE.MAR',
    age = rep(c(20, 80), times = nrow(landscape)),
    biomass = 1
  )

  list(
    landscape = landscape,
    parameters = parameters,
    profile = profile,
    prepared = prepared,
    cohorts = cohorts
  )
}

make_scalar_sequence <- function(values) {
  i <- 0L

  function(n = 1L) {
    if (n != 1L) {
      stop('test sequence supports n=1 only')
    }

    i <<- i + 1L

    if (i > length(values)) {
      stop('test sequence exhausted')
    }

    values[i]
  }
}

test_that('Strict Base Fire profile wires all six process components', {
  fx <- make_basefire_profile_fixture()

  expect_s3_class(
    fx$profile,
    'landisR_fire_basefire_profile'
  )

  expect_s3_class(
    fx$profile$model,
    'landisR_fire_model'
  )

  expect_equal(
    unname(vapply(
      .fire_component_names,
      function(nm) attr(
        fx$profile$model[[nm]],
        'fire_component_class',
        exact = TRUE
      ),
      character(1)
    )),
    unname(.fire_component_names)
  )
})

test_that('Profile preparation builds all six caches once', {
  fx <- make_basefire_profile_fixture()

  expect_s3_class(
    fx$prepared,
    'landisR_fire_basefire_prepared_profile'
  )

  expect_true(all(vapply(
    fx$prepared$model$prepared,
    function(x) !is.null(x),
    logical(1)
  )))
})

test_that('Deterministic complete event wires ignition through mortality', {
  fx <- make_basefire_profile_fixture()

  # Active-site order is 1..9. Only cell 5 passes ignition:
  # cells 1-4: first gate fails because 1 > 0.01 * 10
  # cell 5: ignition + initial-spread gates pass
  # cells 6-9: first gate fails for the same reason
  ignition_draws <- make_scalar_sequence(
    c(
      1, 1, 1, 1,
      0, 0,
      1, 1, 1, 1
    )
  )

  source_id <- fx$landscape$cell_id[
    fx$landscape$row == 2 &
      fx$landscape$col == 2
  ]

  out <- .basefire_run_timestep(
    prepared_profile = fx$prepared,
    current_time = 10,
    cohorts = fx$cohorts,
    rng = list(
      ignition_uniform = ignition_draws,
      exponential_draw = function(lambda) 3,
      next_double = function(n = 1L) rep(0, n),
      generate_uniform = function(n = 1L) rep(0, n),
      shuffle_fun = function(x) x,
      wind_direction = 0L,
      wind_speed = 0
    )
  )

  expect_s3_class(
    out,
    'landisR_fire_basefire_timestep_result'
  )

  expect_equal(out$summary$number_fires, 1L)
  expect_equal(length(out$raw_events), 1L)
  expect_equal(length(out$events), 1L)

  event <- out$events[[1L]]

  expect_s3_class(event, 'landisR_fire_event')
  expect_equal(event$ignition_cell, source_id)
  expect_equal(length(event$cells), 3L)
  expect_equal(event$actual_area, 3)
  expect_equal(event$target_area, 3)
  expect_true(all(event$severity$severity == 2))

  # At severity 2 with tolerance 3, young 20-year cohorts are killed
  # through the first damage class while 80-year cohorts survive.
  expect_equal(
    out$raw_events[[1L]]$cohorts_killed,
    3L
  )

  expect_equal(nrow(out$cohorts), nrow(fx$cohorts) - 3L)

  burned <- out$severity$cell_id[
    out$severity$disturbed
  ]

  expect_equal(length(burned), 3L)
  expect_true(all(
    out$fire_state$cells$last_fire_year[
      match(burned, out$fire_state$cells$cell_id)
    ] == 10
  ))

  expect_equal(nrow(out$event_log), 1L)
  expect_equal(out$event_log$burned_sites, 3L)
  expect_equal(out$event_log$cohorts_killed, 3L)
  expect_equal(out$event_log$severity, 2)

  expect_equal(out$summary$total_sites_burned, 3L)
  expect_equal(unname(out$summary$region_sites['r']), 3L)
})

test_that('Timestep resets event/severity maps but carries fire history', {
  fx <- make_basefire_profile_fixture()

  # No ignition at all.
  out <- .basefire_run_timestep(
    prepared_profile = fx$prepared,
    current_time = 10,
    cohorts = fx$cohorts,
    rng = list(
      ignition_uniform = function(n = 1L) rep(1, n)
    )
  )

  expect_false(any(out$severity$disturbed))
  expect_true(all(out$severity$severity == 0))
  expect_true(all(out$severity$event_owner == 0L))

  initial_last_fire <- fx$prepared$model$prepared$
    ignition$last_fire_year

  expect_equal(
    out$fire_state$cells$last_fire_year,
    initial_last_fire
  )
})

test_that('Dynamic region year requires explicit map reader', {
  fx <- make_basefire_profile_fixture()

  parameters <- fx$parameters
  parameters$dynamic_regions <- data.frame(
    year = 10L,
    map_name = 'dynamic-regions.tif',
    stringsAsFactors = FALSE
  )

  class(parameters) <- 'landisR_fire_basefire_parameters'
  .basefire_validate_parameters(parameters)

  profile <- .basefire_profile(parameters)

  prepared <- .basefire_prepare_profile(
    profile = profile,
    landscape = fx$landscape,
    context = list(
      start_time = 0,
      cell_area = 1,
      species_longevity = c(
        'PICE.MAR' = 100
      )
    )
  )

  expect_error(
    .basefire_run_timestep(
      prepared_profile = prepared,
      current_time = 10,
      cohorts = fx$cohorts,
      rng = list(
        ignition_uniform = function(n = 1L) rep(1, n)
      )
    ),
    'dynamic region map'
  )
})

test_that('Dynamic region reader refreshes region-dependent caches', {
  fx <- make_basefire_profile_fixture()

  parameters <- .basefire_parameters(
    timestep = 10,
    species = data.frame(
      species_code = 'PICE.MAR',
      fire_tolerance = 3
    ),
    regions = data.frame(
      name = c('r1', 'r2'),
      map_code = c(1, 2),
      mean_size = c(3, 3),
      min_size = c(3, 3),
      max_size = c(3, 3),
      ignition_probability = c(0.01, 0.02),
      fire_spread_age = c(20, 40)
    ),
    initial_fire_regions = 'regions.tif',
    dynamic_regions = data.frame(
      year = 10L,
      map_name = 'dynamic-regions.tif'
    ),
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

  profile <- .basefire_profile(parameters)

  prepared <- .basefire_prepare_profile(
    profile = profile,
    landscape = fx$landscape,
    context = list(
      start_time = 0,
      cell_area = 1,
      species_longevity = c(
        'PICE.MAR' = 100
      ),
      dynamic_region_map_reader = function(
          map_name,
          current_time,
          landscape) {
        rep(2L, nrow(landscape))
      }
    )
  )

  out <- .basefire_run_timestep(
    prepared_profile = prepared,
    current_time = 10,
    cohorts = fx$cohorts,
    rng = list(
      ignition_uniform = function(n = 1L) rep(1, n)
    )
  )

  expect_true(all(
    out$prepared_profile$model$prepared$
      ignition$region_index == 2L
  ))

  expect_true(all(
    out$prepared_profile$model$prepared$
      ignition$fire_spread_age == 40
  ))
})

test_that('Base Fire profile remains internal', {
  wd <- normalizePath(
    getwd(),
    winslash = '/',
    mustWork = TRUE
  )

  candidates <- c(
    file.path(wd, 'NAMESPACE'),
    file.path(wd, '..', 'NAMESPACE'),
    file.path(wd, '..', '..', 'NAMESPACE')
  )

  namespace_file <- candidates[file.exists(candidates)][1L]

  expect_false(is.na(namespace_file))

  namespace <- readLines(
    namespace_file,
    warn = FALSE
  )

  expect_false(any(grepl(
    'basefire',
    namespace,
    ignore.case = TRUE
  )))
})
