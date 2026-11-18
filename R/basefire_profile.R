.basefire_model <- function() {
  fire_model(
    ignition = .basefire_ignition,
    size = .basefire_size,
    spread = .basefire_spread,
    burn_probability = .basefire_burn_probability,
    severity = .basefire_severity,
    mortality = .basefire_mortality,
    name = 'basefire-strict',
    metadata = list(
      profile = 'basefire',
      compatibility = 'strict',
      public = TRUE
    )
  )
}

.basefire_profile <- function(parameters) {
  .basefire_validate_parameters(parameters)

  x <- list(
    parameters = parameters,
    timestep = as.integer(parameters$timestep),
    model = .basefire_model()
  )

  class(x) <- 'landisR_fire_basefire_profile'
  .basefire_validate_profile(x)
  x
}

.basefire_validate_profile <- function(x) {
  if (!inherits(x, 'landisR_fire_basefire_profile')) {
    .fire_stop(
      'x must inherit from landisR_fire_basefire_profile.'
    )
  }

  .basefire_validate_parameters(x$parameters)
  validate_fire_model(x$model)

  if (!identical(x$model$name, 'basefire-strict')) {
    .fire_stop('Invalid Base Fire strict model name.')
  }

  if (!identical(
    unname(vapply(
      .fire_component_names,
      function(nm) attr(
        x$model[[nm]],
        'fire_component_class',
        exact = TRUE
      ),
      character(1)
    )),
    unname(.fire_component_names)
  )) {
    .fire_stop('Base Fire profile component wiring is invalid.')
  }

  invisible(TRUE)
}

.basefire_prepare_profile <- function(
    profile,
    landscape,
    fire_state = NULL,
    context = list(),
    topology = NULL,
    site_topology = NULL) {

  .basefire_validate_profile(profile)
  .assert_named_list(context, 'context')

  if (!is.data.frame(landscape)) {
    .fire_stop('Base Fire profile landscape must be a data.frame.')
  }

  parameter_cache <- .basefire_prepare_parameter_cache(
    profile$parameters
  )

  context$basefire_cache <- parameter_cache

  # Strict GetNeighbors() needs the complete geometric landscape, including
  # inactive sites. Always refresh this pointer to the current landscape.
  context$basefire_full_landscape <- landscape

  prepare_context <- context
  prepare_context$basefire_site_topology <- site_topology
  prepared_model <- fire_prepare(
    model = profile$model,
    landscape = landscape,
    fire_state = fire_state,
    context = prepare_context,
    topology = topology
  )
  # The spread cache owns the site topology, including in saved checkpoints.
  prepared_model$prepare_context$basefire_site_topology <- NULL

  x <- list(
    profile = profile,
    model = prepared_model,
    parameters = profile$parameters,
    parameter_cache = parameter_cache,
    landscape = landscape,
    fire_state = fire_state,
    context = context,
    topology = prepared_model$topology
  )

  class(x) <- 'landisR_fire_basefire_prepared_profile'
  .basefire_validate_prepared_profile(x)
  x
}

.basefire_validate_prepared_profile <- function(x) {
  if (!inherits(
    x,
    'landisR_fire_basefire_prepared_profile'
  )) {
    .fire_stop(
      'x must inherit from ',
      'landisR_fire_basefire_prepared_profile.'
    )
  }

  .basefire_validate_profile(x$profile)
  validate_fire_prepared_model(x$model)

  if (!inherits(
    x$parameter_cache,
    'landisR_fire_basefire_parameter_cache'
  )) {
    .fire_stop('Invalid Base Fire shared parameter cache.')
  }

  for (nm in .fire_component_names) {
    if (is.null(x$model$prepared[[nm]])) {
      .fire_stop(
        'Base Fire prepared profile is missing cache: ',
        nm,
        '.'
      )
    }
  }

  invisible(TRUE)
}

.basefire_validate_cohort_table_legacy_270B <- function(
    cohorts,
    topology) {

  if (!is.data.frame(cohorts)) {
    .fire_stop('cohorts must be a data.frame.')
  }

  .basefire_require_columns(
    cohorts,
    c('cell_id', 'species_code', 'age'),
    'cohorts'
  )

  if (nrow(cohorts) == 0L) {
    return(invisible(TRUE))
  }

  cell_id <- .basefire_integer_like(
    cohorts$cell_id,
    'cohorts$cell_id'
  )

  species_code <- .basefire_nonempty_character(
    as.character(cohorts$species_code),
    'cohorts$species_code'
  )

  if (!is.numeric(cohorts$age) ||
      anyNA(cohorts$age) ||
      any(!is.finite(cohorts$age)) ||
      any(cohorts$age < 0)) {
    .fire_stop(
      'cohorts$age must contain finite values >= 0.'
    )
  }

  unknown_cell <- setdiff(
    unique(cell_id),
    topology$cell_id
  )

  if (length(unknown_cell)) {
    .fire_stop(
      'cohorts reference inactive/unknown cell_id values: ',
      paste(unknown_cell, collapse = ', '),
      '.'
    )
  }

  invisible(TRUE)
}

.basefire_normalize_rng <- function(rng = list()) {
  .assert_named_list(rng, 'rng')

  defaults <- list(
    ignition_uniform = stats::runif,
    generate_uniform = stats::runif,
    next_double = stats::runif,
    exponential_draw = .basefire_default_exponential_draw,
    shuffle_fun = .basefire_default_shuffle,
    wind_direction = NULL,
    wind_speed = NULL
  )

  unknown <- setdiff(names(rng), names(defaults))

  if (length(unknown)) {
    .fire_stop(
      'Unknown Base Fire RNG controls: ',
      paste(unknown, collapse = ', '),
      '.'
    )
  }

  out <- defaults
  out[names(rng)] <- rng

  for (nm in c(
    'ignition_uniform',
    'generate_uniform',
    'next_double',
    'exponential_draw',
    'shuffle_fun'
  )) {
    .assert_function(out[[nm]], paste0('rng$', nm))
  }

  if (!is.null(out$wind_direction)) {
    if (!is.numeric(out$wind_direction) ||
        length(out$wind_direction) != 1L ||
        is.na(out$wind_direction) ||
        !is.finite(out$wind_direction) ||
        out$wind_direction < 0 ||
        abs(out$wind_direction - round(out$wind_direction)) >
          sqrt(.Machine$double.eps)) {
      .fire_stop(
        'rng$wind_direction must be NULL or one integer >= 0.'
      )
    }

    out$wind_direction <- as.integer(
      round(out$wind_direction)
    )
  }

  if (!is.null(out$wind_speed)) {
    .basefire_wind_probabilities(out$wind_speed)
  }

  out
}

.basefire_apply_mortality_component <- function(
    cohorts,
    severity,
    mortality_fun,
    mortality_prepared) {

  .assert_function(mortality_fun, 'mortality_fun')

  if (!nrow(cohorts)) {
    return(list(
      cohorts_killed = 0L,
      kill = logical(),
      removal_fraction = numeric()
    ))
  }

  kill <- logical(nrow(cohorts))
  removal_fraction <- numeric(nrow(cohorts))

  for (i in seq_len(nrow(cohorts))) {
    fraction <- mortality_fun(
      severity = severity,
      species_code = as.character(
        cohorts$species_code[i]
      ),
      cohort_age = cohorts$age[i],
      prepared = mortality_prepared,
      diagnostics = FALSE
    )

    if (!is.numeric(fraction) ||
        length(fraction) != 1L ||
        is.na(fraction) ||
        !is.finite(fraction) ||
        fraction < 0 ||
        fraction > 1) {
      .fire_stop(
        'Base Fire mortality component returned ',
        'an invalid removal fraction.'
      )
    }

    removal_fraction[i] <- fraction
    kill[i] <- fraction >= 1
  }

  list(
    cohorts_killed = as.integer(sum(kill)),
    kill = kill,
    removal_fraction = removal_fraction
  )
}

.basefire_spread_to_fire_event <- function(
    spread_result,
    prepared_profile,
    current_time,
    event_id) {

  burned <- spread_result$cells[
    spread_result$cells$burned,
    ,
    drop = FALSE
  ]

  if (!nrow(burned)) {
    return(NULL)
  }

  initial_region <- prepared_profile$model$prepared$
    spread$region_index[
      spread_result$initiation_cell_index
    ]

  cell_area <- prepared_profile$model$prepared$size$cell_area

  target_area <-
    spread_result$size_targets[initial_region] *
    cell_area

  actual_area <-
    spread_result$total_sites_in_event *
    cell_area

  fire_event(
    year = current_time,
    ignition_cell = spread_result$initiation_cell_id,
    cells = burned$cell_id,
    severity = data.frame(
      cell_id = burned$cell_id,
      severity = burned$severity
    ),
    target_area = target_area,
    actual_area = actual_area,
    fire_region = prepared_profile$parameter_cache$
      region_name[initial_region],
    metadata = list(
      profile = 'basefire-strict',
      event_id = as.integer(event_id),
      sites_checked = spread_result$num_site_checked,
      cohorts_killed = spread_result$cohorts_killed,
      average_event_severity = spread_result$severity,
      size_targets = spread_result$size_targets,
      sites_in_event = spread_result$sites_in_event,
      wind_direction_raw =
        spread_result$wind_direction_raw,
      wind_direction = spread_result$wind_direction,
      wind_speed = spread_result$wind_speed
    )
  )
}

.basefire_event_log_row <- function(
    spread_result,
    prepared_profile,
    current_time,
    event_id) {

  # PlugIn.LogEvent() writes only when FireEvent.Severity > 0.
  if (!(spread_result$severity > 0)) {
    return(NULL)
  }

  landscape <- prepared_profile$landscape
  ignition_id <- spread_result$initiation_cell_id
  row_index <- match(
    ignition_id,
    landscape$cell_id
  )

  data.frame(
    event_id = as.integer(event_id),
    time = as.integer(current_time),
    row = as.integer(landscape$row[row_index]),
    column = as.integer(landscape$col[row_index]),
    sites_checked = as.integer(
      spread_result$num_site_checked
    ),
    cohorts_killed = as.integer(
      spread_result$cohorts_killed
    ),
    severity = as.numeric(spread_result$severity),
    burned_sites = as.integer(
      spread_result$total_sites_in_event
    ),
    stringsAsFactors = FALSE
  )
}

.basefire_dynamic_region_codes <- function(
    map_object,
    landscape) {

  if (is.data.frame(map_object)) {
    .basefire_require_columns(
      map_object,
      c('cell_id', 'fire_region_map_code'),
      'dynamic fire-region map'
    )

    idx <- match(
      landscape$cell_id,
      map_object$cell_id
    )

    if (anyNA(idx)) {
      .fire_stop(
        'Dynamic fire-region map is missing landscape cells.'
      )
    }

    codes <- map_object$fire_region_map_code[idx]
  } else {
    codes <- map_object
  }

  codes <- .basefire_integer_like(
    codes,
    'dynamic fire-region map codes'
  )

  if (length(codes) != nrow(landscape)) {
    .fire_stop(
      'Dynamic fire-region map must provide one code ',
      'per landscape row.'
    )
  }

  if (any(codes < 0L | codes > 65535L)) {
    .fire_stop(
      'Dynamic fire-region map codes must be 0..65535.'
    )
  }

  codes
}

.basefire_apply_dynamic_region_update <- function(
    prepared_profile,
    current_time) {

  years <- prepared_profile$parameter_cache$dynamic_year
  hit <- which(years == current_time)

  if (!length(hit)) {
    return(prepared_profile)
  }

  if (length(hit) != 1L) {
    .fire_stop(
      'Multiple dynamic Base Fire region maps for one year.'
    )
  }

  reader <- prepared_profile$context$
    dynamic_region_map_reader

  if (is.null(reader) || !is.function(reader)) {
    .fire_stop(
      'A Base Fire dynamic region map is scheduled for year ',
      current_time,
      ' but context$dynamic_region_map_reader ',
      'is not available.'
    )
  }

  map_name <- prepared_profile$parameter_cache$
    dynamic_map_name[hit]

  timer <- .fire_timer("carte_dynamique", current_time)
  map_object <- reader(
    map_name = map_name,
    current_time = current_time,
    landscape = prepared_profile$landscape
  )

  timer("lecture")
  codes <- .basefire_dynamic_region_codes(
    map_object,
    prepared_profile$landscape
  )

  landscape <- prepared_profile$landscape
  landscape$fire_region_map_code <- codes

  # Dynamic maps replace the fire-region map. Avoid a stale explicit
  # fire_region_index column overriding the new map codes.
  if ('fire_region_index' %in% names(landscape)) {
    landscape$fire_region_index <- NULL
  }

  timer("alignement_codes")
  # Region codes change independently of cell order, coordinates and mask.
  updated <- .basefire_prepare_profile(
    profile = prepared_profile$profile,
    landscape = landscape,
    fire_state = prepared_profile$fire_state,
    context = prepared_profile$context,
    topology = prepared_profile$topology,
    site_topology = prepared_profile$model$prepared$spread$site_topology
  )
  timer("reconstruction_modele")
  updated
}

.basefire_run_timestep_legacy_270C <- function(
    prepared_profile,
    current_time,
    cohorts,
    rng = list(), report = NULL) {

  .basefire_validate_prepared_profile(
    prepared_profile
  )

  if (!is.numeric(current_time) ||
      length(current_time) != 1L ||
      is.na(current_time) ||
      !is.finite(current_time)) {
    .fire_stop(
      'current_time must be one finite numeric value.'
    )
  }

  current_time <- as.integer(current_time)
  if (!is.null(report)) report(0.13, "Fire regions & cohort index")

  prepared_profile <-
    .basefire_apply_dynamic_region_update(
      prepared_profile,
      current_time
    )

  topology <- prepared_profile$topology

  .basefire_validate_cohort_table(
    cohorts,
    topology
  )

  rng <- .basefire_normalize_rng(rng)

  cohort_env <- new.env(parent = emptyenv())
  cohort_env$table <- cohorts
  rownames(cohort_env$table) <- NULL

  n <- length(topology$cell_id)

  last_fire_year <- prepared_profile$model$prepared$
    ignition$last_fire_year

  last_severity <- rep(NA_real_, n)

  if (!is.null(prepared_profile$fire_state)) {
    state_idx <- match(
      topology$cell_id,
      prepared_profile$fire_state$cells$cell_id
    )

    if (anyNA(state_idx)) {
      .fire_stop(
        'Prepared Base Fire state is missing active cells.'
      )
    }

    last_severity <-
      prepared_profile$fire_state$cells$
        last_severity[state_idx]
  }

  runtime <- new.env(parent = emptyenv())
  runtime$last_fire_year <- as.numeric(last_fire_year)
  runtime$last_severity <- as.numeric(last_severity)
  runtime$event_occupied <- rep(FALSE, n)
  runtime$event_owner <- integer(n)
  runtime$severity <- numeric(n)
  runtime$disturbed <- rep(FALSE, n)

  last_wind_year <- prepared_profile$model$prepared$
    severity$last_wind_year

  severity_fun <- function(cell_index, current_time) {
    prepared_profile$model$severity(
      cell_index = cell_index,
      current_time = current_time,
      prepared = prepared_profile$model$prepared$severity,
      last_fire_year = runtime$last_fire_year,
      last_wind_year = last_wind_year,
      diagnostics = FALSE
    )
  }

  species_codes <- prepared_profile$model$prepared$mortality$species_code
  species_tally <- matrix(0, nrow = length(species_codes), ncol = 2L,
                          dimnames = list(species_codes, c("cohorts", "density")))
  damage_fun <- function(
      cell_index,
      severity,
      current_time) {

    cell_id <- topology$cell_id[cell_index]

    cohort_rows <- which(
      cohort_env$table$cell_id == cell_id
    )

    if (!length(cohort_rows)) {
      return(0L)
    }

    local <- cohort_env$table[
      cohort_rows,
      ,
      drop = FALSE
    ]

    mortality <- .basefire_apply_mortality_component(
      cohorts = local,
      severity = severity,
      mortality_fun = prepared_profile$model$mortality,
      mortality_prepared =
        prepared_profile$model$prepared$mortality
    )

    if (any(mortality$kill)) {
      global_kill <- cohort_rows[mortality$kill]
      killed <- local[mortality$kill, , drop = FALSE]
      cohort_env$killed_density <- cohort_env$killed_density + sum(killed$biomass)
      tally <- rowsum(cbind(cohorts = rep.int(1L, nrow(killed)), density = killed$biomass),
                      killed$species_code)
      positions <- match(rownames(tally), species_codes)
      cohort_env$killed_species[positions, ] <-
        cohort_env$killed_species[positions, , drop = FALSE] + tally
      cohort_env$table <- cohort_env$table[
        -global_kill,
        ,
        drop = FALSE
      ]
      rownames(cohort_env$table) <- NULL
    }

    mortality$cohorts_killed
  }

  raw_events <- list()
  events <- list()
  log_rows <- list()
  killed_density <- numeric()
  killed_species <- list()

  summary_region_sites <- integer(
    prepared_profile$parameter_cache$n_regions
  )

  summary_total_sites <- 0L
  summary_event_count <- 0L

  stride <- if (!is.null(report)) max(1L, ceiling(n / 20))
  for (cell_index in seq_len(n)) {
    if (!is.null(report) && (cell_index - 1L) %% stride == 0)
      report(0.23 + 0.48 * (cell_index - 1L) / n, "Ignition & spread")
    # Keep ignition TSF state synchronized with events that have already
    # occurred earlier in this same timestep.
    ignition_prepared <-
      prepared_profile$model$prepared$ignition

    ignition_prepared$last_fire_year <-
      runtime$last_fire_year

    initiation <- prepared_profile$model$ignition(
      cell_index = cell_index,
      current_time = current_time,
      timestep = prepared_profile$profile$timestep,
      severity_fun = severity_fun,
      prepared = ignition_prepared,
      uniform = rng$ignition_uniform,
      diagnostics = TRUE
    )

    if (!isTRUE(initiation$initiated)) {
      next
    }

    event_id <- length(raw_events) + 1L
    cohort_env$killed_density <- 0
    cohort_env$killed_species <- species_tally

    spread_result <- .basefire_run_spread_event(
      initiation_cell_index = cell_index,
      current_time = current_time,
      spread_prepared =
        prepared_profile$model$prepared$spread,
      size_prepared =
        prepared_profile$model$prepared$size,
      burn_prepared =
        prepared_profile$model$prepared$burn_probability,
      severity_fun = severity_fun,
      damage_fun = damage_fun,
      size_fun = prepared_profile$model$size,
      spread_fun = prepared_profile$model$spread,
      burn_probability_fun =
        prepared_profile$model$burn_probability,
      exponential_draw = rng$exponential_draw,
      generate_uniform = rng$generate_uniform,
      next_double = rng$next_double,
      shuffle_fun = rng$shuffle_fun,
      event_occupied = runtime$event_occupied,
      last_fire_year = runtime$last_fire_year,
      wind_direction = rng$wind_direction,
      wind_speed = rng$wind_speed
    )

    runtime$last_fire_year <-
      spread_result$last_fire_year

    runtime$event_occupied <-
      spread_result$event_occupied

    member <- which(spread_result$event_member)

    if (length(member)) {
      # Strict cross-event behavior: a newly initiated event can overwrite
      # ownership of its initiation site. Other occupied neighbors cannot join.
      runtime$event_owner[member] <- event_id
    }

    burned <- spread_result$cells[
      spread_result$cells$burned,
      ,
      drop = FALSE
    ]

    if (nrow(burned)) {
      burned_index <- burned$cell_index

      runtime$severity[burned_index] <-
        burned$severity

      runtime$disturbed[burned_index] <- TRUE
      runtime$last_severity[burned_index] <-
        burned$severity
    }

    raw_events[[event_id]] <- spread_result
    summary_event_count <- summary_event_count + 1L

    standard_event <- .basefire_spread_to_fire_event(
      spread_result = spread_result,
      prepared_profile = prepared_profile,
      current_time = current_time,
      event_id = event_id
    )

    if (!is.null(standard_event)) {
      events[[length(events) + 1L]] <- standard_event
    }

    log_row <- .basefire_event_log_row(
      spread_result = spread_result,
      prepared_profile = prepared_profile,
      current_time = current_time,
      event_id = event_id
    )

    if (!is.null(log_row)) {
      log_rows[[length(log_rows) + 1L]] <- log_row
      killed_density[length(log_rows)] <- cohort_env$killed_density
      tally <- cohort_env$killed_species
      killed_species[[length(log_rows)]] <- tally[tally[, "cohorts"] > 0, , drop = FALSE]

      summary_total_sites <-
        summary_total_sites +
        spread_result$total_sites_in_event

      summary_region_sites <-
        summary_region_sites +
        spread_result$sites_in_event
    }
  }

  if (!is.null(report)) report(0.73, "Building fire logs & history")
  if (length(log_rows)) {
    event_log <- do.call(rbind, log_rows)
    rownames(event_log) <- NULL
  } else {
    event_log <- data.frame(
      event_id = integer(),
      time = integer(),
      row = integer(),
      column = integer(),
      sites_checked = integer(),
      cohorts_killed = integer(),
      severity = numeric(),
      burned_sites = integer(),
      stringsAsFactors = FALSE
    )
  }

  names(summary_region_sites) <-
    prepared_profile$parameter_cache$region_name

  summary <- list(
    time = current_time,
    total_sites_burned = as.integer(
      summary_total_sites
    ),
    number_fires = as.integer(
      summary_event_count
    ),
    region_sites = summary_region_sites
  )

  state <- fire_state(
    cell_id = topology$cell_id,
    last_fire_year = runtime$last_fire_year,
    last_severity = runtime$last_severity
  )

  if (!is.null(last_wind_year)) {
    state$cells$last_wind_year <- last_wind_year
  }

  # Carry the updated state into the prepared profile so the same object can
  # be used for the next timestep without reinitializing fire history.
  prepared_profile$fire_state <- state

  prepared_profile$model$prepared$ignition$
    last_fire_year <- runtime$last_fire_year

  prepared_profile$model$prepared$burn_probability$
    last_fire_year <- runtime$last_fire_year

  prepared_profile$model$prepared$severity$
    last_fire_year <- runtime$last_fire_year

  out <- list(
    time = current_time,
    events = events,
    raw_events = raw_events,
    event_log = event_log,
    killed_density = killed_density,
    killed_species = killed_species,
    summary = summary,
    cohorts = cohort_env$table,
    fire_state = state,
    severity = data.frame(
      cell_id = topology$cell_id,
      severity = runtime$severity,
      disturbed = runtime$disturbed,
      event_owner = runtime$event_owner
    ),
    prepared_profile = prepared_profile
  )

  class(out) <- 'landisR_fire_basefire_timestep_result'
  out
}
