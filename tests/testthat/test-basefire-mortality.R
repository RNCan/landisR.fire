make_mortality_cache <- function(
    tolerance = c('PICE.MAR' = 3L, 'PINU.BAN' = 2L),
    longevity = c('PICE.MAR' = 100, 'PINU.BAN' = 200),
    max_age = c(0.25, 0.50, 1.00),
    damage_diff = c(-2L, 1L, 4L)) {

  species <- names(tolerance)

  x <- list(
    species_code = species,
    species_fire_tolerance = as.integer(unname(tolerance)),
    species_longevity = as.numeric(unname(longevity[species])),
    species_index = setNames(seq_along(species), species),
    damage_max_age = max_age,
    damage_severity_tolerance_diff = damage_diff
  )

  class(x) <- 'landisR_fire_basefire_mortality_cache'
  x
}

test_that('Severity 5 kills all cohorts unconditionally', {
  prepared <- make_mortality_cache()

  out <- .basefire_cohort_mortality(
    severity = 5,
    species_code = 'UNKNOWN',
    cohort_age = 100000,
    prepared = prepared,
    diagnostics = TRUE
  )

  expect_true(out$killed)
  expect_equal(out$removal_fraction, 1)
  expect_equal(out$reason, 'severity5')
  expect_true(is.na(out$fire_tolerance))
  expect_true(is.na(out$longevity))
})

test_that('Damage table uses severity minus species fire tolerance', {
  prepared <- make_mortality_cache()

  # severity 2 - tolerance 3 = -1.
  # First damage class requires >= -2 and age <= 25% longevity.
  young <- .basefire_cohort_mortality(
    severity = 2,
    species_code = 'PICE.MAR',
    cohort_age = 20,
    prepared = prepared,
    diagnostics = TRUE
  )

  old <- .basefire_cohort_mortality(
    severity = 2,
    species_code = 'PICE.MAR',
    cohort_age = 30,
    prepared = prepared,
    diagnostics = TRUE
  )

  expect_true(young$killed)
  expect_equal(young$matched_damage_index, 1L)
  expect_equal(young$severity_tolerance_difference, -1L)

  expect_false(old$killed)
  expect_equal(old$removal_fraction, 0)
})

test_that('Damage MaxAge comparison is inclusive', {
  prepared <- make_mortality_cache()

  at_boundary <- .basefire_cohort_mortality(
    severity = 2,
    species_code = 'PICE.MAR',
    cohort_age = 25,
    prepared = prepared,
    diagnostics = TRUE
  )

  expect_true(at_boundary$killed)
  expect_equal(at_boundary$age_as_percent, .basefire_float32(0.25))
})

test_that('Damage severity-difference comparison is inclusive', {
  prepared <- make_mortality_cache()

  # severity 4 - tolerance 3 = 1, exactly the second damage threshold.
  out <- .basefire_cohort_mortality(
    severity = 4,
    species_code = 'PICE.MAR',
    cohort_age = 50,
    prepared = prepared,
    diagnostics = TRUE
  )

  expect_true(out$killed)
  expect_equal(out$severity_tolerance_difference, 1L)
  expect_equal(out$matched_damage_index, 2L)
})

test_that('Cohort older than longevity survives non-severity-5 fire when unmatched', {
  prepared <- make_mortality_cache()

  out <- .basefire_cohort_mortality(
    severity = 4,
    species_code = 'PICE.MAR',
    cohort_age = 120,
    prepared = prepared,
    diagnostics = TRUE
  )

  expect_gt(out$age_as_percent, 1)
  expect_false(out$killed)
})

test_that('Executed forward damage-table order is preserved', {
  prepared <- make_mortality_cache(
    max_age = c(0.25, 0.50, 1.00),
    damage_diff = c(-2L, 1L, 4L)
  )

  # Young cohort at severity 4 satisfies both class 1 and class 2.
  # Current C# foreach sees class 1 first because parser appends in increasing
  # order and PlugIn passes FireDamages unchanged.
  out <- .basefire_cohort_mortality(
    severity = 4,
    species_code = 'PICE.MAR',
    cohort_age = 20,
    prepared = prepared,
    diagnostics = TRUE
  )

  expect_true(out$killed)
  expect_equal(out$matched_damage_index, 1L)
})

test_that('Mortality returns exactly 1 or 0 like disturbance callback', {
  prepared <- make_mortality_cache()

  killed <- .basefire_cohort_mortality(
    severity = 2,
    species_code = 'PICE.MAR',
    cohort_age = 20,
    prepared = prepared
  )

  survives <- .basefire_cohort_mortality(
    severity = 2,
    species_code = 'PICE.MAR',
    cohort_age = 80,
    prepared = prepared
  )

  expect_identical(killed, 1.0)
  expect_identical(survives, 0.0)
})

test_that('Site mortality returns kill mask and cohort count', {
  prepared <- make_mortality_cache()

  cohorts <- data.frame(
    species_code = c(
      'PICE.MAR',
      'PICE.MAR',
      'PINU.BAN',
      'PINU.BAN'
    ),
    age = c(20, 80, 40, 190),
    biomass = c(10, 20, 30, 40)
  )

  out <- .basefire_apply_mortality(
    cohorts = cohorts,
    severity = 2,
    prepared = prepared
  )

  expect_equal(out$cohorts_killed, 2L)
  expect_equal(out$kill, c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(out$removal_fraction, c(1, 0, 1, 0))
  expect_equal(nrow(out$survivors), 2L)
  expect_equal(nrow(out$killed), 2L)
})

test_that('Float32 helper reproduces single-precision boundary behavior', {
  x <- 1 / 3
  x32 <- .basefire_float32(x)

  expect_true(is.numeric(x32))
  expect_length(x32, 1L)
  expect_false(identical(x32, x))
  expect_equal(
    .basefire_float32(x32),
    x32
  )
})

test_that('Mortality preparation aligns fire tolerance and longevity', {
  parameters <- .basefire_parameters(
    timestep = 10,
    species = data.frame(
      species_code = c('PICE.MAR', 'PINU.BAN'),
      fire_tolerance = c(3, 2)
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
    damage = data.frame(
      max_age = c(0.25, 0.5, 1),
      severity_tolerance_diff = c(-2, 1, 4)
    ),
    map_names_template = 'severity.tif',
    log_file = 'events.csv',
    summary_log_file = 'summary.csv'
  )

  cache <- .basefire_prepare_parameter_cache(parameters)

  prepared <- .basefire_prepare_mortality(
    landscape = NULL,
    fire_state = NULL,
    topology = NULL,
    context = list(
      basefire_cache = cache,
      species_longevity = c(
        'PINU.BAN' = 200,
        'PICE.MAR' = 100
      )
    ),
    model = NULL
  )

  expect_s3_class(
    prepared,
    'landisR_fire_basefire_mortality_cache'
  )

  expect_equal(
    prepared$species_code,
    c('PICE.MAR', 'PINU.BAN')
  )

  expect_equal(
    prepared$species_fire_tolerance,
    c(3L, 2L)
  )

  expect_equal(
    prepared$species_longevity,
    c(100, 200)
  )

  expect_equal(
    prepared$damage_max_age,
    c(0.25, 0.5, 1)
  )

  expect_equal(
    prepared$damage_severity_tolerance_diff,
    c(-2L, 1L, 4L)
  )
})

test_that('Mortality preparation also accepts species table', {
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
    damage = data.frame(
      max_age = 1,
      severity_tolerance_diff = 4
    ),
    map_names_template = 'severity.tif',
    log_file = 'events.csv',
    summary_log_file = 'summary.csv'
  )

  cache <- .basefire_prepare_parameter_cache(parameters)

  prepared <- .basefire_prepare_mortality(
    landscape = NULL,
    fire_state = NULL,
    topology = NULL,
    context = list(
      basefire_cache = cache,
      species = data.frame(
        species_code = 'PICE.MAR',
        longevity = 150
      )
    ),
    model = NULL
  )

  expect_equal(prepared$species_longevity, 150)
})

test_that('Base Fire mortality is a real internal mortality component', {
  expect_true(is.function(.basefire_mortality))

  expect_equal(
    attr(.basefire_mortality, 'fire_component_class'),
    'mortality'
  )

  expect_true(is.function(
    attr(.basefire_mortality, 'fire_prepare')
  ))
})
