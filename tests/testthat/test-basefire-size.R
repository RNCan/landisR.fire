make_basefire_size_cache <- function(
    mean_size = 100,
    min_size = 10,
    max_size = 1000,
    cell_area = 1) {

  x <- list(
    region_index = 1L,
    basefire_index = 0L,
    region_name = 'region-1',
    mean_size = mean_size,
    min_size = min_size,
    max_size = max_size,
    cell_area = cell_area
  )

  class(x) <- 'landisR_fire_basefire_size_cache'
  x
}

make_exponential_sequence <- function(values) {
  i <- 0L

  function(lambda) {
    i <<- i + 1L

    if (i > length(values)) {
      stop('exponential sequence exhausted')
    }

    values[i]
  }
}

test_that('MeanSize <= 0 returns zero without consuming a draw', {
  prepared <- make_basefire_size_cache(
    mean_size = 0,
    min_size = 50,
    max_size = 1000,
    cell_area = 10
  )

  draw_calls <- 0L
  draw <- function(lambda) {
    draw_calls <<- draw_calls + 1L
    123
  }

  out <- .basefire_compute_size(
    region_index = 1,
    prepared = prepared,
    exponential_draw = draw,
    diagnostics = TRUE
  )

  expect_equal(out$size_sites, 0L)
  expect_equal(out$draws, 0L)
  expect_equal(draw_calls, 0L)
  expect_equal(out$bound, 'mean_nonpositive')
})

test_that('ComputeSize sets lambda to inverse MeanSize', {
  prepared <- make_basefire_size_cache(
    mean_size = 200,
    min_size = 10,
    max_size = 1000,
    cell_area = 1
  )

  observed_lambda <- NA_real_

  draw <- function(lambda) {
    observed_lambda <<- lambda
    200
  }

  out <- .basefire_compute_size(
    1,
    prepared,
    exponential_draw = draw,
    diagnostics = TRUE
  )

  expect_equal(observed_lambda, 1 / 200)
  expect_equal(out$lambda, 1 / 200)
  expect_equal(out$draws, 1L)
})

test_that('Generated fire size below MinSize is replaced by MinSize', {
  prepared <- make_basefire_size_cache(
    mean_size = 100,
    min_size = 25,
    max_size = 1000,
    cell_area = 10
  )

  out <- .basefire_compute_size(
    1,
    prepared,
    exponential_draw = function(lambda) 5,
    diagnostics = TRUE
  )

  expect_equal(out$size_generated, 5)
  expect_equal(out$bounded_size, 25)
  expect_equal(out$size_sites, 2L)
  expect_equal(out$bound, 'min')
})

test_that('Generated fire size above MaxSize is replaced by MaxSize', {
  prepared <- make_basefire_size_cache(
    mean_size = 100,
    min_size = 25,
    max_size = 250,
    cell_area = 100
  )

  out <- .basefire_compute_size(
    1,
    prepared,
    exponential_draw = function(lambda) 500,
    diagnostics = TRUE
  )

  expect_equal(out$bounded_size, 250)
  expect_equal(out$size_sites, 2L)
  expect_equal(out$bound, 'max')
})

test_that('ComputeSize truncates positive cell count toward zero', {
  prepared <- make_basefire_size_cache(
    mean_size = 100,
    min_size = 0,
    max_size = 1000,
    cell_area = 10
  )

  expect_equal(
    .basefire_compute_size(
      1,
      prepared,
      exponential_draw = function(lambda) 29.999
    ),
    2L
  )

  expect_equal(
    .basefire_compute_size(
      1,
      prepared,
      exponential_draw = function(lambda) 30
    ),
    3L
  )
})

test_that('Base Fire does not force a positive fire to at least one cell', {
  prepared <- make_basefire_size_cache(
    mean_size = 100,
    min_size = 1,
    max_size = 1000,
    cell_area = 100
  )

  out <- .basefire_compute_size(
    1,
    prepared,
    exponential_draw = function(lambda) 5,
    diagnostics = TRUE
  )

  expect_equal(out$bounded_size, 5)
  expect_equal(out$size_sites, 0L)
})

test_that('Strict less-than and greater-than bounds are preserved', {
  prepared <- make_basefire_size_cache(
    mean_size = 100,
    min_size = 25,
    max_size = 250,
    cell_area = 10
  )

  at_min <- .basefire_compute_size(
    1,
    prepared,
    exponential_draw = function(lambda) 25,
    diagnostics = TRUE
  )

  at_max <- .basefire_compute_size(
    1,
    prepared,
    exponential_draw = function(lambda) 250,
    diagnostics = TRUE
  )

  expect_equal(at_min$bound, 'none')
  expect_equal(at_max$bound, 'none')
  expect_equal(at_min$size_sites, 2L)
  expect_equal(at_max$size_sites, 25L)
})

test_that('Unordered Min/Mean/Max values are not silently normalized', {
  prepared <- make_basefire_size_cache(
    mean_size = 20,
    min_size = 50,
    max_size = 10,
    cell_area = 1
  )

  # This mirrors the literal C# branch order. A generated value of 30 first
  # satisfies generated < MinSize and therefore becomes 50, even though MaxSize
  # is 10.
  out <- .basefire_compute_size(
    1,
    prepared,
    exponential_draw = function(lambda) 30,
    diagnostics = TRUE
  )

  expect_equal(out$bounded_size, 50)
  expect_equal(out$bound, 'min')
  expect_equal(out$size_sites, 50L)
})

test_that('Size preparation extracts compact region vectors and cell area', {
  parameters <- .basefire_parameters(
    timestep = 10,
    species = data.frame(
      species_code = 'PICE.MAR',
      fire_tolerance = 3
    ),
    regions = data.frame(
      name = c('north', 'south'),
      map_code = c(1, 2),
      mean_size = c(100, 200),
      min_size = c(10, 20),
      max_size = c(1000, 2000),
      ignition_probability = c(0.01, 0.02),
      fire_spread_age = c(20, 30)
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

  prepared <- .basefire_prepare_size(
    landscape = NULL,
    fire_state = NULL,
    topology = NULL,
    context = list(
      basefire_cache = cache,
      cell_area = 6.25
    ),
    model = NULL
  )

  expect_s3_class(prepared, 'landisR_fire_basefire_size_cache')
  expect_equal(prepared$mean_size, c(100, 200))
  expect_equal(prepared$min_size, c(10, 20))
  expect_equal(prepared$max_size, c(1000, 2000))
  expect_equal(prepared$cell_area, 6.25)
})

test_that('Base Fire size is now a real internal size component', {
  expect_true(is.function(.basefire_size))
  expect_equal(
    attr(.basefire_size, 'fire_component_class'),
    'size'
  )
  expect_true(is.function(
    attr(.basefire_size, 'fire_prepare')
  ))
})
