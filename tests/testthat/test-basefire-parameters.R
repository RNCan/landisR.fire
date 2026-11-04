test_that('Base Fire regions reproduce reference validation rules', {
  regions <- data.frame(
    name = c('north', 'south'),
    map_code = c(1, 65535),
    mean_size = c(100, 20),
    min_size = c(5, 50),
    max_size = c(1000, 10),
    ignition_probability = c(0.01, 1),
    fire_spread_age = c(20, 0)
  )

  out <- .basefire_normalize_regions(regions)

  expect_equal(out$region_index, 1:2)
  expect_equal(out$basefire_index, 0:1)
  expect_equal(out$map_code, c(1L, 65535L))

  # Deliberately accepted for strict compatibility: the C# parser does not
  # impose MinSize <= MeanSize <= MaxSize.
  expect_equal(out$min_size[2], 50)
  expect_equal(out$max_size[2], 10)

  bad_name <- regions
  bad_name$name[2] <- 'north'
  expect_error(.basefire_normalize_regions(bad_name), 'unique')

  bad_code <- regions
  bad_code$map_code[2] <- 1
  expect_error(.basefire_normalize_regions(bad_code), 'unique')

  bad_prob <- regions
  bad_prob$ignition_probability[1] <- 1.1
  expect_error(.basefire_normalize_regions(bad_prob), 'between 0 and 1')
})

test_that('Fuel and wind curve rules remain distinct', {
  regions <- .basefire_normalize_regions(data.frame(
    name = c('north', 'south'),
    map_code = c(1, 2),
    mean_size = c(100, 100),
    min_size = c(5, 5),
    max_size = c(1000, 1000),
    ignition_probability = c(0.01, 0.02),
    fire_spread_age = c(20, 30)
  ))

  fuel <- data.frame(
    region_name = 'north',
    severity1 = -1,
    severity2 = 5,
    severity3 = 10,
    severity4 = 20,
    severity5 = 30
  )
  expect_silent(.basefire_normalize_curve(fuel, regions, 'fuel'))

  fuel_equal <- fuel
  fuel_equal$severity1 <- 5
  fuel_equal$severity2 <- 5
  expect_error(
    .basefire_normalize_curve(fuel_equal, regions, 'fuel'),
    'must be -1 or >'
  )

  wind_equal <- data.frame(
    region_name = 'north',
    severity1 = 0,
    severity2 = 0,
    severity3 = 0,
    severity4 = 0,
    severity5 = 0
  )
  expect_silent(.basefire_normalize_curve(wind_equal, regions, 'wind'))
})

test_that('Prepared Base Fire parameter cache is compact and direct', {
  p <- .basefire_parameters(
    timestep = 10,
    species = data.frame(
      species_code = c('PICE.MAR', 'PINU.BAN'),
      fire_tolerance = c(3, 4)
    ),
    regions = data.frame(
      name = c('north', 'south'),
      map_code = c(2, 65535),
      mean_size = c(100, 200),
      min_size = c(5, 10),
      max_size = c(1000, 2000),
      ignition_probability = c(0.01, 0.02),
      fire_spread_age = c(20, 30)
    ),
    initial_fire_regions = 'regions.tif',
    dynamic_regions = data.frame(
      year = c(50, 100),
      map_name = c('regions50.tif', 'regions100.tif')
    ),
    fuel_curves = data.frame(
      region_name = 'north',
      severity1 = -1,
      severity2 = 5,
      severity3 = 10,
      severity4 = 20,
      severity5 = 30
    ),
    wind_curves = data.frame(
      region_name = c('north', 'south'),
      severity1 = c(0, -1),
      severity2 = c(0, -1),
      severity3 = c(0, -1),
      severity4 = c(0, -1),
      severity5 = c(0, -1)
    ),
    damage = data.frame(
      max_age = c(0.25, 0.5, 1),
      severity_tolerance_diff = c(-2, 1, 4)
    ),
    map_names_template = 'severity-{timestep}.tif',
    log_file = 'fire-events.csv',
    summary_log_file = 'fire-summary.csv'
  )

  cache <- .basefire_prepare_parameter_cache(p)

  expect_s3_class(cache, 'landisR_fire_basefire_parameter_cache')
  expect_length(cache$map_code_to_region, 65536)
  expect_equal(cache$map_code_to_region[2 + 1L], 1L)
  expect_equal(cache$map_code_to_region[65535 + 1L], 2L)

  expect_equal(unname(cache$fuel_curve['north', ]), c(-1L, 5L, 10L, 20L, 30L))
  expect_equal(unname(cache$fuel_curve['south', ]), rep(0L, 5))
  expect_equal(unname(cache$species_fire_tolerance['PINU.BAN']), 4L)
})
