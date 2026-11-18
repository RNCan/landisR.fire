# 215-internal-handoff-v1
# Base Fire strict <-> landisR adapter.
# Low-level helpers remain internal; public wrappers are defined below.

.basefire_require_landisR <- function() {
  if (!requireNamespace('landisR', quietly = TRUE)) {
    .fire_stop(
      'Base Fire landisR coupling requires landisR >= 0.1.0.9008.'
    )
  }

  version <- utils::packageVersion('landisR')

  if (version < package_version('0.1.0.9008')) {
    .fire_stop(
      'Base Fire landisR coupling requires landisR >= 0.1.0.9008; ',
      'installed version is ',
      as.character(version),
      '.'
    )
  }

  invisible(version)
}

.basefire_landisR_csr <- function(simulation) {
  .basefire_require_landisR()
  if (!"landis_csr" %in% getNamespaceExports("landisR")) return(NULL)
  view <- landisR::landis_csr(simulation)
  if (identical(view$states$format, "landisR52/state-csr-v3")) view
}

.basefire_landisR_state <- function(simulation) {
  view <- .basefire_landisR_csr(simulation)
  if (!is.null(view)) {
    return(basefire_expand_csr(view$states, view$cell_state, view$cell_id, view$species))
  }
  cohorts <- landisR::landis_state(
    simulation,
    representation = 'cohorts'
  )

  required <- c('cell_id', 'species', 'age', 'biomass')
  missing <- setdiff(required, names(cohorts))

  if (length(missing)) {
    .fire_stop(
      'landis_state(cohorts) is missing required columns: ',
      paste(missing, collapse = ', '),
      '.'
    )
  }

  out <- data.frame(
    cell_id = as.integer(cohorts$cell_id),
    species_code = as.character(cohorts$species),
    age = as.numeric(cohorts$age),
    biomass = as.numeric(cohorts$biomass),
    stringsAsFactors = FALSE
  )

  rownames(out) <- NULL
  out
}

.basefire_align_fire_regions <- function(value, cells) {
  n <- nrow(cells)

  if (is.data.frame(value)) {
    .basefire_require_columns(
      value,
      c('cell_id', 'fire_region_map_code'),
      'fire_region_map_code'
    )

    idx <- match(cells$cell_id, value$cell_id)

    if (anyNA(idx)) {
      .fire_stop(
        'fire_region_map_code data.frame must contain every landscape cell.'
      )
    }

    value <- value$fire_region_map_code[idx]
  } else if (length(value) == 1L) {
    value <- rep(value, n)
  } else if (!is.null(names(value)) && all(nzchar(names(value)))) {
    idx <- match(as.character(cells$cell_id), names(value))
    if (anyNA(idx)) {
      .fire_stop(
        'Named fire_region_map_code vector must contain every landscape cell.'
      )
    }
    value <- value[idx]
  }

  value <- .basefire_integer_like(value, 'fire_region_map_code')

  if (length(value) != n) {
    .fire_stop(
      'fire_region_map_code must have length 1 or one value per landscape cell.'
    )
  }

  value
}

.basefire_landisR_landscape <- function(
    simulation,
    fire_region_map_code) {

  .basefire_require_landisR()

  if (!inherits(simulation, 'landisR_simulation')) {
    .fire_stop('simulation must inherit from landisR_simulation.')
  }

  canonical <- simulation$project$canonical
  cells <- canonical$landscape$cells

  required <- c('cell_id', 'row', 'col', 'active')
  missing <- setdiff(required, names(cells))

  if (length(missing)) {
    .fire_stop(
      'landisR canonical landscape is missing columns: ',
      paste(missing, collapse = ', '),
      '.'
    )
  }

  max_age <- numeric(nrow(cells))
  view <- .basefire_landisR_csr(simulation)
  if (is.null(view)) {
    state <- .basefire_landisR_state(simulation)
    if (nrow(state)) {
      split_age <- split(state$age, state$cell_id)
      ids <- as.integer(names(split_age))
      max_age[match(ids, cells$cell_id)] <- vapply(split_age, max, numeric(1))
    }
  } else {
    pos <- match(view$cell_id, cells$cell_id)
    if (anyNA(pos))
      .fire_stop("Base Fire CSR cell identifiers do not match the landscape.")
    max_age[pos] <- basefire_csr_max_age(view$states, view$cell_state)
  }

  data.frame(
    cell_id = as.integer(cells$cell_id),
    row = as.integer(cells$row),
    col = as.integer(cells$col),
    active = as.logical(cells$active),
    fire_region_map_code = .basefire_align_fire_regions(
      fire_region_map_code,
      cells
    ),
    max_cohort_age = as.numeric(max_age),
    stringsAsFactors = FALSE
  )
}

.basefire_prepare_landisR <- function(
    simulation,
    parameters,
    fire_region_map_code,
    cell_area,
    fire_state = NULL,
    dynamic_region_map_reader = NULL,
    time_of_last_wind = NULL,
    shared_geometry = FALSE) {

  timer <- .fire_timer("initialisation", simulation$year)
  .basefire_require_landisR()
  .basefire_validate_parameters(parameters)

  if (!inherits(simulation, 'landisR_simulation')) {
    .fire_stop('simulation must inherit from landisR_simulation.')
  }

  if (!is.numeric(cell_area) ||
      length(cell_area) != 1L ||
      is.na(cell_area) ||
      !is.finite(cell_area) ||
      cell_area <= 0) {
    .fire_stop('cell_area must be one finite numeric value > 0.')
  }

  canonical <- simulation$project$canonical

  if (!identical(
    as.integer(parameters$timestep),
    as.integer(canonical$timestep)
  )) {
    .fire_stop(
      'Base Fire timestep (',
      parameters$timestep,
      ') must equal landisR succession timestep (',
      canonical$timestep,
      ') for the current coupled adapter.'
    )
  }

  species <- canonical$species

  if (!all(c('species', 'longevity') %in% names(species))) {
    .fire_stop(
      'landisR species table must contain species and longevity.'
    )
  }

  missing_species <- setdiff(
    parameters$species$species_code,
    as.character(species$species)
  )

  if (length(missing_species)) {
    .fire_stop(
      'Base Fire species absent from landisR: ',
      paste(missing_species, collapse = ', '),
      '.'
    )
  }

  longevity <- as.numeric(species$longevity)
  names(longevity) <- as.character(species$species)

  grid <- canonical$landscape
  raster_reader <- NULL
  if (!is.null(parameters$region_rasters)) {
    raster_reader <- .basefire_raster_map_reader(
      parameters$region_rasters, grid$nrow, grid$ncol
    )
  }
  if (is.null(fire_region_map_code)) {
    fire_region_map_code <- parameters$initial_region_codes
    if (is.null(fire_region_map_code)) {
      if (is.null(raster_reader))
        .fire_stop('Supply fire_region_map_code or build parameters with basefire_build().')
      fire_region_map_code <- raster_reader(parameters$initial_fire_regions, landscape = grid$cells)
    }
  }
  if (is.null(dynamic_region_map_reader)) {
    if (length(parameters$dynamic_region_maps))
      dynamic_region_map_reader <- .basefire_table_map_reader(parameters$dynamic_region_maps)
    else if (nrow(parameters$dynamic_regions))
      dynamic_region_map_reader <- raster_reader
  }

  timer("validation")
  landscape <- .basefire_landisR_landscape(
    simulation = simulation,
    fire_region_map_code = fire_region_map_code
  )

  timer("vue_paysage")
  context <- list(
    start_time = as.integer(simulation$year),
    cell_area = as.numeric(cell_area),
    species_longevity = longevity
  )

  if (!is.null(dynamic_region_map_reader)) {
    .assert_function(
      dynamic_region_map_reader,
      'dynamic_region_map_reader'
    )
    context$dynamic_region_map_reader <- dynamic_region_map_reader
  }

  if (!is.null(time_of_last_wind)) {
    context$time_of_last_wind <- time_of_last_wind
  }

  profile <- .basefire_profile(parameters)

  timer("parametres_et_contexte")
  geometry <- if (shared_geometry) .basefire_shared_geometry(landscape)
  prepared <- .basefire_prepare_profile(
    profile = profile,
    landscape = landscape,
    fire_state = fire_state,
    context = context,
    topology = geometry$active,
    site_topology = geometry$sites
  )
  timer("preparation_modele")
  prepared
}

.basefire_landisR_update <- function(
    fire_result,
    burned_cells, preserve_cohort_order = FALSE) {

  .basefire_require_landisR()

  burned_cells <- unique(as.integer(burned_cells))

  if (!length(burned_cells)) {
    return(NULL)
  }

  survivors <- fire_result$cohorts[
    fire_result$cohorts$cell_id %in% burned_cells,
    ,
    drop = FALSE
  ]

  if (nrow(survivors)) {
    replacement <- data.frame(
      cell_id = as.integer(survivors$cell_id),
      species = as.character(survivors$species_code),
      age = as.numeric(survivors$age),
      biomass = as.numeric(survivors$biomass),
      stringsAsFactors = FALSE
    )
  } else {
    replacement <- NULL
  }

  nonempty <- if (is.null(replacement)) {
    integer()
  } else {
    unique(as.integer(replacement$cell_id))
  }

  empty <- setdiff(burned_cells, nonempty)

  landisR::landis_update(
    cohorts = replacement,
    clear_cells = empty,
    preserve_cohort_order = preserve_cohort_order,
    metadata = list(
      source = 'landisR.fire',
      profile = 'basefire-strict',
      process = 'mortality'
    )
  )
}

.basefire_step_landisR <- function(
    simulation,
    prepared_profile,
    rng = list(),
    stochastic_serotiny = TRUE,
    metadata = list(),
    year = simulation$year,
    use_landis_rng = FALSE,
    report_progress = getOption("landisR.progress.report")) {

  timer <- .fire_timer("pas", year)
  .basefire_require_landisR()
  .basefire_validate_prepared_profile(prepared_profile)

  if (!inherits(simulation, 'landisR_simulation')) {
    .fire_stop('simulation must inherit from landisR_simulation.')
  }

  if (!is.logical(stochastic_serotiny) ||
      length(stochastic_serotiny) != 1L ||
      is.na(stochastic_serotiny)) {
    .fire_stop('stochastic_serotiny must be TRUE or FALSE.')
  }

  .assert_named_list(metadata, 'metadata')

  if (!is.numeric(year) || length(year) != 1L || is.na(year) ||
      !year %in% c(simulation$year, simulation$year + simulation$project$canonical$timestep)) {
    .fire_stop('Fire year must be the current or next succession year.')
  }

  timer("validation")
  if (!is.null(report_progress)) report_progress(0.06, "Reading cohorts")
  cohorts <- .basefire_landisR_state(simulation)
  timer("expansion_cohortes")

  fire_result <- .basefire_run_timestep(
    prepared_profile = prepared_profile,
    current_time = as.integer(year),
    cohorts = cohorts,
    rng = rng,
    landis_rng = if (use_landis_rng) simulation$landis_rng,
    report = report_progress
  )
  timer("moteur")
  if (use_landis_rng) simulation$landis_rng <- fire_result$landis_rng

  burned <- fire_result$severity[
    fire_result$severity$disturbed &
      fire_result$severity$severity > 0,
    c('cell_id', 'severity'),
    drop = FALSE
  ]

  timer("cellules_brulees")
  if (!nrow(burned)) {
    out <- list(
      simulation = simulation,
      fire_result = fire_result,
      prepared_profile = fire_result$prepared_profile,
      update = NULL,
      burned_cells = integer()
    )
    class(out) <- 'landisR_fire_basefire_landisR_result'
    return(out)
  }

  if (!is.null(report_progress)) report_progress(0.84, "Preparing cohort mortality")
  update <- .basefire_landisR_update(
    fire_result = fire_result,
    burned_cells = burned$cell_id,
    preserve_cohort_order = use_landis_rng
  )

  timer("preparation_mortalite")
  if (!is.null(report_progress)) report_progress(0.89, "Applying mortality & serotiny")
  simulation <- landisR::landis_fire(
    simulation = simulation,
    cells = as.integer(burned$cell_id),
    update = update,
    stochastic = stochastic_serotiny,
    severity = as.numeric(burned$severity),
    year = year,
    metadata = c(
      list(
        source = 'landisR.fire',
        profile = 'basefire-strict',
        basefire_events = length(fire_result$events)
      ),
      metadata
    )
  )

  timer("application_mortalite_regeneration")
  out <- list(
    simulation = simulation,
    fire_result = fire_result,
    prepared_profile = fire_result$prepared_profile,
    update = update,
    burned_cells = as.integer(burned$cell_id)
  )

  class(out) <- 'landisR_fire_basefire_landisR_result'
  out
}


# 216B-public-api-v1
# Public high-level Base Fire strict API.

basefire_read <- function(
    path,
    known_species = NULL,
    base_dir = dirname(path),
    project = NULL,
    scenario_directory = dirname(path),
    raster_template = NULL) {

  parameters <- .basefire_read_parameters(
    path = path,
    known_species = known_species,
    base_dir = base_dir,
    project = project
  )
  .assert_scalar_character(scenario_directory, 'scenario_directory')
  directory <- normalizePath(scenario_directory, winslash = '/', mustWork = TRUE)
  template <- if (is.null(raster_template)) parameters$initial_fire_regions else raster_template
  .assert_scalar_character(template, 'raster_template')
  parameters$region_rasters <- list(
    directory = directory,
    template = .basefire_resolve_path(chartr('\\', '/', template), directory)
  )
  parameters
}

basefire_prepare <- function(
    simulation,
    parameters,
    fire_region_map_code = NULL,
    cell_area,
    fire_state = NULL,
    dynamic_region_map_reader = NULL,
    time_of_last_wind = NULL) {

  .basefire_prepare_landisR(
    simulation = simulation,
    parameters = parameters,
    fire_region_map_code = fire_region_map_code,
    cell_area = cell_area,
    fire_state = fire_state,
    dynamic_region_map_reader = dynamic_region_map_reader,
    time_of_last_wind = time_of_last_wind
  )
}

.basefire_step_legacy_270B <- function(
    simulation,
    prepared,
    rng = list(),
    stochastic_serotiny = TRUE,
    metadata = list()) {

  .basefire_step_landisR(
    simulation = simulation,
    prepared_profile = prepared,
    rng = rng,
    stochastic_serotiny = stochastic_serotiny,
    metadata = metadata
  )
}
