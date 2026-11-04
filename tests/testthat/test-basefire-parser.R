test_that('Current Base Fire text format parses to normalized parameters', {
  td <- tempfile('basefire-')
  dir.create(td)

  species_file <- file.path(td, 'species.csv')
  writeLines(
    c(
      'SpeciesCode,FireTolerance',
      'PICE.MAR,3',
      'PINU.BAN,4'
    ),
    species_file
  )

  input_file <- file.path(td, 'base-fire.txt')
  writeLines(
    c(
      'LandisData "Base Fire"',
      'Timestep 10',
      'Species_CSV_File species.csv',
      '',
      '>> Region MapCode Mean Min Max Ignition FireSpreadAge',
      'north 1 100 5 1000 0.01 20',
      'south 2 200 10 2000 0.02 30',
      'InitialFireRegionsMap regions.tif',
      '',
      'DynamicFireRegionTable',
      '50 regions50.tif',
      '100 regions100.tif',
      '',
      'FuelCurveTable',
      'north -1 5 10 20 30',
      '',
      'WindCurveTable',
      # Current C# reads file order as S5 S4 S3 S2 S1.
      'north 40 30 20 10 0',
      'south -1 -1 -1 -1 -1',
      '',
      'FireDamageTable',
      '25% -2',
      '50% 1',
      '100% 4',
      'MapNames severity-{timestep}.tif',
      'LogFile fire-events.csv',
      'SummaryLogFile fire-summary.csv'
    ),
    input_file
  )

  p <- .basefire_read_parameters(
    input_file,
    known_species = c('PICE.MAR', 'PINU.BAN', 'ABIE.BAL')
  )

  expect_s3_class(p, 'landisR_fire_basefire_parameters')
  expect_equal(p$timestep, 10L)
  expect_equal(p$regions$name, c('north', 'south'))
  expect_equal(p$dynamic_regions$year, c(50L, 100L))
  expect_equal(p$damage$max_age, c(0.25, 0.5, 1))
  expect_equal(p$wind_curves$severity1[1], 0L)
  expect_equal(p$wind_curves$severity5[1], 40L)
  expect_equal(attr(p, 'landis_data'), 'Base Fire')

  cache <- .basefire_prepare_parameter_cache(p)
  expect_equal(cache$region_by_name[['south']], 2L)
})

test_that('Base Fire parser rejects unknown supplied species', {
  td <- tempfile('basefire-species-')
  dir.create(td)

  species_file <- file.path(td, 'species.csv')
  writeLines(
    c(
      'SpeciesCode,FireTolerance',
      'UNKNOWN,3'
    ),
    species_file
  )

  input_file <- file.path(td, 'base-fire.txt')
  writeLines(
    c(
      'LandisData "Base Fire"',
      'Timestep 10',
      'Species_CSV_File species.csv',
      'north 1 100 5 1000 0.01 20',
      'InitialFireRegionsMap regions.tif',
      'FuelCurveTable',
      'WindCurveTable',
      'FireDamageTable',
      '100% 4',
      'MapNames severity.tif',
      'LogFile fire-events.csv',
      'SummaryLogFile fire-summary.csv'
    ),
    input_file
  )

  expect_error(
    .basefire_read_parameters(input_file, known_species = 'PICE.MAR'),
    'not present'
  )
})
