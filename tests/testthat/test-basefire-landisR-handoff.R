# 215-real-landisR-handoff-v1

# Exact fixture copied from landisR's validated public round-trip test.
make_serotiny_project <- function(type = "strict", dispersal = "local", source_age = 30L, ncol = 5L, include_traits = TRUE) {
    species <- data.frame(species = "PINU.BAN", longevity = 120L, maturity = 15L, shade_tolerance = 1L, fire_tolerance = 3L,
        effective_seed_distance = 100, max_seed_distance = 250, vegetative_probability = 0, min_sprout_age = 0L, max_sprout_age = 0L,
        post_fire_regeneration = "serotiny", stringsAsFactors = FALSE)
    if (isTRUE(include_traits)) {
        species$serotiny_type <- type
        species$serotiny_dispersal <- dispersal
    }
    species_parameters <- data.frame(species = "PINU.BAN", leaf_longevity = 3, woody_decay = 0.05, mortality_shape = 15,
        growth_shape = 0.5, leaf_lignin = 0.2, stringsAsFactors = FALSE)
    ecoregions <- data.frame(map_code = 1L, ecoregion = "eco1", active = TRUE, description = "test", stringsAsFactors = FALSE)
    ecoregion_parameters <- data.frame(ecoregion = "eco1", aet = 600, stringsAsFactors = FALSE)
    dynamic <- data.frame(year = 0L, ecoregion = "eco1", species = "PINU.BAN", sep = 1, max_anpp = 100, max_biomass = 5000,
        stringsAsFactors = FALSE)
    initial_cohorts <- data.frame(map_code = 1L, species = "PINU.BAN", age = source_age, biomass = 100, stringsAsFactors = FALSE)
    shade <- data.frame(shade = 1:5, eco1 = c(0, 0.2, 0.4, 0.6, 0.8))
    light <- data.frame(shade_tolerance = 1:5, shade_0 = 1, shade_1 = 1, shade_2 = 1, shade_3 = 1, shade_4 = 1, shade_5 = 1)
    fire_reductions <- data.frame(severity = 1L, wood_reduction = 0, litter_reduction = 0.5)
    landscape <- data.frame(cell_id = seq_len(ncol), row = rep(1L, ncol), col = seq_len(ncol), ecoregion_map_code = rep(1L,
        ncol), ecoregion = rep("eco1", ncol), active = rep(TRUE, ncol), community = c(1L, rep(2L, ncol - 1L)), stringsAsFactors = FALSE)
    landisR::landis_build(landscape = landscape, species = species, species_parameters = species_parameters, ecoregions = ecoregions,
        ecoregion_parameters = ecoregion_parameters, dynamic = dynamic, initial_cohorts = initial_cohorts, shade = shade,
        light = light, fire_reductions = fire_reductions, duration = 10L, cell_length = 100, timestep = 10L, seeding_algorithm = "WardSeedDispersal",
        spinup_cohorts = FALSE, spinup_mortality_fraction = 0, initial_mode = "explicit_biomass", community_codes = c(1L,
            2L), name = "serotiny-test")
}

make_215_scalar_sequence <- function(values) {
  i <- 0L
  function(n = 1L) {
    if (!identical(as.integer(n), 1L)) {
      stop('215 deterministic sequence supports n=1 only')
    }
    i <<- i + 1L
    if (i > length(values)) stop('215 deterministic sequence exhausted')
    values[i]
  }
}

test_that('Base Fire strict handoff runs on a real landisR simulation', {
  skip_if_not_installed(
    'landisR',
    minimum_version = '0.1.0.9003'
  )

  project <- make_serotiny_project(
    type = 'strict',
    dispersal = 'local',
    ncol = 5L,
    include_traits = TRUE
  )

  sim0 <- landisR::landis_initialize(
    project,
    seed = 4357L
  )

  pre <- landisR::landis_state(
    sim0,
    representation = 'cohorts'
  )

  species_row <- project$canonical$species[
    project$canonical$species$species == 'PINU.BAN',
    ,
    drop = FALSE
  ]

  expect_equal(nrow(species_row), 1L)

  mature <- pre[
    pre$species == 'PINU.BAN' &
      pre$age >= species_row$maturity[1L],
    ,
    drop = FALSE
  ]

  expect_gt(nrow(mature), 0L)

  target_cell <- as.integer(mature$cell_id[1L])
  active_ids <- as.integer(sim0$ctx$ids)
  target_index <- match(target_cell, active_ids)

  expect_false(is.na(target_index))

  parameters <- .basefire_parameters(
    timestep = project$canonical$timestep,
    species = data.frame(
      species_code = 'PINU.BAN',
      fire_tolerance = as.integer(species_row$fire_tolerance[1L]),
      stringsAsFactors = FALSE
    ),
    regions = data.frame(
      name = 'r1',
      map_code = 1L,
      mean_size = 1,
      min_size = 1,
      max_size = 1,
      ignition_probability = 0.01,
      fire_spread_age = 20L,
      stringsAsFactors = FALSE
    ),
    initial_fire_regions = 'internal-test-map',
    dynamic_regions = NULL,
    fuel_curves = data.frame(
      region_name = 'r1',
      severity1 = 0L,
      severity2 = 1L,
      severity3 = 2L,
      severity4 = 3L,
      severity5 = 4L,
      stringsAsFactors = FALSE
    ),
    wind_curves = NULL,
    damage = data.frame(
      max_age = 1,
      severity_tolerance_diff = 4L,
      stringsAsFactors = FALSE
    ),
    map_names_template = 'severity-{timestep}.tif',
    log_file = 'events.csv',
    summary_log_file = 'summary.csv'
  )

  prepared <- basefire_prepare(
    simulation = sim0,
    parameters = parameters,
    fire_region_map_code = 1L,
    cell_area = 1
  )

  n_active <- length(active_ids)
  ignition_values <- c(
    if (target_index > 1L) rep(1, target_index - 1L) else numeric(),
    0,
    0,
    if (target_index < n_active) rep(1, n_active - target_index) else numeric()
  )

  out <- basefire_step(
    simulation = sim0,
    prepared = prepared,
    rng = list(
      ignition_uniform = make_215_scalar_sequence(ignition_values),
      exponential_draw = function(lambda) 1,
      next_double = function(n = 1L) rep(1, n),
      generate_uniform = function(n = 1L) rep(1, n),
      shuffle_fun = function(x) x,
      wind_direction = 0L,
      wind_speed = 0
    ),
    stochastic_serotiny = FALSE,
    metadata = list(test = '215')
  )

  expect_s3_class(
    out,
    'landisR_fire_basefire_landisR_result'
  )

  expect_identical(out$burned_cells, target_cell)
  expect_equal(out$fire_result$summary$number_fires, 1L)
  expect_equal(out$fire_result$summary$total_sites_burned, 1L)

  burned_row <- out$fire_result$severity[
    out$fire_result$severity$cell_id == target_cell,
    ,
    drop = FALSE
  ]

  expect_true(burned_row$disturbed)
  expect_equal(burned_row$severity, 5)

  # Severity 5 removes all pre-fire cohorts in the target cell.
  expect_false(any(
    out$fire_result$cohorts$cell_id == target_cell
  ))

  expect_true(target_cell %in% out$update$clear_cells)

  # landis_fire() must capture pre-fire serotinous sources BEFORE applying
  # the mortality update, then recruit from those sources.
  post <- landisR::landis_state(
    out$simulation,
    representation = 'cohorts'
  )

  recruits <- post[
    post$cell_id == target_cell &
      post$species == 'PINU.BAN' &
      post$age == 1,
    ,
    drop = FALSE
  ]

  expect_gt(nrow(recruits), 0L)

  # Untargeted cells must be unchanged by the mortality handoff.
  pre_other <- pre[pre$cell_id != target_cell, , drop = FALSE]
  post_other <- post[post$cell_id != target_cell, , drop = FALSE]

  ord_pre <- order(pre_other$cell_id, pre_other$species, pre_other$age)
  ord_post <- order(post_other$cell_id, post_other$species, post_other$age)

  rownames(pre_other) <- NULL
  rownames(post_other) <- NULL

  expect_identical(
    pre_other[ord_pre, , drop = FALSE],
    post_other[ord_post, , drop = FALSE]
  )

  expect_equal(
    out$prepared_profile$fire_state$cells$last_fire_year[
      match(target_cell, out$prepared_profile$fire_state$cells$cell_id)
    ],
    sim0$year
  )
})

test_that('Base Fire reads scenario-relative rasters and applies a scheduled map', {
  skip_if_not_installed('terra')
  skip_if_not_installed('landisR')
  skip_if_not_installed('landisR.succession.biomass')

  td <- tempfile('basefire-rasters-')
  scenario <- file.path(td, 'replicate')
  dir.create(scenario, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  project <- make_serotiny_project(ncol = 5L)
  grid <- terra::rast(nrow = 1, ncol = 5, xmin = 0, xmax = 500,
                     ymin = 0, ymax = 100, crs = 'EPSG:32618')
  initial <- c(1L, 2L, 1L, 1L, 2L)
  dynamic <- c(2L, 1L, 2L, 2L, 1L)
  terra::writeRaster(terra::setValues(grid, initial), file.path(td, 'initial.tif'))
  dynamic_file <- file.path(td, 'dynamic.tif')
  terra::writeRaster(terra::setValues(grid, dynamic), dynamic_file)
  input <- file.path(td, 'base-fire.txt')
  writeLines(c(
    'LandisData "Base Fire"', 'Timestep 10',
    'north 1 1 1 1 0 20', 'south 2 1 1 1 0 20',
    'InitialFireRegionsMap ../initial.tif',
    'DynamicFireRegionTable', '10 ../dynamic.tif',
    'FuelCurveTable', 'WindCurveTable', 'FireDamageTable', '100% 4',
    'MapNames severity.tif', 'LogFile fire.csv', 'SummaryLogFile summary.csv'
  ), input)
  parameters <- basefire_read(input, project = project, scenario_directory = scenario)
  sim <- landisR::landis_initialize(project, seed = 4357L)
  auto <- landisR::landis_sim(sim, steps = 0L, extensions = list(
    fire = basefire_extension(parameters)
  ))
  expect_identical(auto$extension_states$fire$prepared$landscape$fire_region_map_code, initial)

  manual <- landisR::landis_sim(sim, steps = 1L, extensions = list(
    fire = basefire_extension(parameters, fire_region_map_code = initial,
                             dynamic_region_map_reader = function(...) dynamic)
  ))
  resumed <- landisR::landis_sim(unserialize(serialize(auto, NULL)), steps = 1L)
  expect_identical(resumed$extension_states$fire$prepared$landscape$fire_region_map_code, dynamic)
  expect_identical(resumed$events$fire, manual$events$fire)
  expect_identical(landisR::landis_csr(resumed$simulation), landisR::landis_csr(manual$simulation))

  terra::writeRaster(terra::shift(terra::setValues(grid, dynamic), dx = 100),
                     dynamic_file, overwrite = TRUE)
  expect_error(landisR::landis_sim(auto, steps = 1L), 'not aligned to the reference grid')
})
