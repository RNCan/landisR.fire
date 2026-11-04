.basefire_compute_init_spread_prob <- function(
    time_since_last_fire,
    fire_spread_age) {

  if (!is.numeric(time_since_last_fire) ||
      !is.numeric(fire_spread_age) ||
      length(time_since_last_fire) != length(fire_spread_age) &&
        length(time_since_last_fire) != 1L &&
        length(fire_spread_age) != 1L) {
    .fire_stop(
      'time_since_last_fire and fire_spread_age must be numeric and ',
      'vector-compatible.'
    )
  }

  # Intentionally mirrors FireEvent.ComputeFireInitSpreadProb exactly:
  #
  #   1 - exp(timeSinceLastFire * (-1 / FireSpreadAge))
  #
  # No special-case clamp is added for FireSpreadAge == 0. In IEEE floating
  # point arithmetic this yields:
  #   TSF > 0, age = 0  -> 1
  #   TSF = 0, age = 0  -> NaN
  #
  # The latter fails a <= comparison, matching the C# behavior.
  1.0 - exp(time_since_last_fire * (-1.0 / fire_spread_age))
}

.basefire_initial_last_fire_year <- function(max_cohort_age, start_time) {
  if (!is.numeric(max_cohort_age) ||
      any(!is.finite(max_cohort_age)) ||
      any(max_cohort_age < 0)) {
    .fire_stop('max_cohort_age must contain finite values >= 0.')
  }

  if (!is.numeric(start_time) ||
      length(start_time) != 1L ||
      is.na(start_time) ||
      !is.finite(start_time)) {
    .fire_stop('start_time must be one finite numeric value.')
  }

  start_time - max_cohort_age
}

.basefire_resolve_region_index <- function(
    landscape,
    topology,
    cache) {

  if (!inherits(cache, 'landisR_fire_basefire_parameter_cache')) {
    .fire_stop(
      'context$basefire_cache must inherit from ',
      'landisR_fire_basefire_parameter_cache.'
    )
  }

  if (!is.data.frame(landscape)) {
    .fire_stop('Base Fire ignition preparation requires a data.frame landscape.')
  }

  if (is.null(topology)) {
    if (!all(c('cell_id', 'row', 'col') %in% names(landscape))) {
      .fire_stop(
        'A topology or landscape columns cell_id, row and col are required.'
      )
    }
    topology <- fire_grid_topology(landscape)
  }

  validate_fire_topology(topology)

  row_index <- match(topology$cell_id, landscape$cell_id)
  if (anyNA(row_index)) {
    .fire_stop('Topology cells are not all present in landscape.')
  }

  if ('fire_region_index' %in% names(landscape)) {
    region_index <- as.integer(landscape$fire_region_index[row_index])
  } else if ('fire_region_map_code' %in% names(landscape)) {
    map_code <- as.integer(landscape$fire_region_map_code[row_index])

    if (anyNA(map_code) || any(map_code < 0L | map_code > 65535L)) {
      .fire_stop(
        'fire_region_map_code must contain integers between 0 and 65535.'
      )
    }

    region_index <- cache$map_code_to_region[map_code + 1L]
  } else {
    .fire_stop(
      'landscape must contain fire_region_index or fire_region_map_code.'
    )
  }

  if (anyNA(region_index) ||
      any(region_index < 1L | region_index > cache$n_regions)) {
    .fire_stop('Landscape contains undefined Base Fire region values.')
  }

  region_index
}

.basefire_prepare_ignition <- function(
    landscape,
    fire_state,
    topology,
    context,
    model) {

  cache <- context$basefire_cache

  if (is.null(cache)) {
    .fire_stop(
      'Base Fire ignition preparation requires context$basefire_cache.'
    )
  }

  if (is.null(topology)) {
    topology <- fire_grid_topology(landscape)
  }

  region_index <- .basefire_resolve_region_index(
    landscape = landscape,
    topology = topology,
    cache = cache
  )

  n <- length(topology$cell_id)

  if (is.null(fire_state)) {
    last_fire_year <- rep(NA_real_, n)
  } else {
    validate_fire_state(fire_state)

    state_index <- match(
      topology$cell_id,
      fire_state$cells$cell_id
    )

    if (anyNA(state_index)) {
      .fire_stop(
        'fire_state is missing at least one active topology cell.'
      )
    }

    last_fire_year <- as.numeric(
      fire_state$cells$last_fire_year[state_index]
    )
  }

  missing_last_fire <- is.na(last_fire_year)

  if (any(missing_last_fire)) {
    if (!('max_cohort_age' %in% names(landscape))) {
      .fire_stop(
        'Base Fire requires last_fire_year or landscape$max_cohort_age ',
        'for initialization.'
      )
    }

    start_time <- context$start_time
    if (is.null(start_time)) {
      .fire_stop(
        'Base Fire initialization from max_cohort_age requires ',
        'context$start_time.'
      )
    }

    landscape_index <- match(topology$cell_id, landscape$cell_id)
    max_age <- landscape$max_cohort_age[landscape_index]

    initialized <- .basefire_initial_last_fire_year(
      max_cohort_age = max_age,
      start_time = start_time
    )

    last_fire_year[missing_last_fire] <- initialized[missing_last_fire]
  }

  x <- list(
    cell_id = topology$cell_id,
    region_index = region_index,
    ignition_probability = cache$ignition_probability[region_index],
    fire_spread_age = cache$fire_spread_age[region_index],
    last_fire_year = last_fire_year
  )

  class(x) <- 'landisR_fire_basefire_ignition_cache'
  x
}

.basefire_validate_ignition_cache <- function(x) {
  if (!inherits(x, 'landisR_fire_basefire_ignition_cache')) {
    .fire_stop(
      'prepared must inherit from landisR_fire_basefire_ignition_cache.'
    )
  }

  n <- length(x$cell_id)

  required <- c(
    'region_index',
    'ignition_probability',
    'fire_spread_age',
    'last_fire_year'
  )

  for (nm in required) {
    if (length(x[[nm]]) != n) {
      .fire_stop('Invalid Base Fire ignition cache field: ', nm, '.')
    }
  }

  invisible(TRUE)
}

.basefire_uniform_one <- function(uniform) {
  .assert_function(uniform, 'uniform')

  value <- uniform(1L)

  if (!is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value < 0 ||
      value > 1) {
    .fire_stop('uniform(1) must return one finite value between 0 and 1.')
  }

  as.numeric(value)
}

.basefire_try_initiate <- function(
    cell_index,
    current_time,
    timestep,
    prepared,
    severity_fun,
    uniform = stats::runif,
    diagnostics = FALSE) {

  .basefire_validate_ignition_cache(prepared)
  .assert_function(severity_fun, 'severity_fun')
  .assert_function(uniform, 'uniform')

  cell_index <- as.integer(cell_index)

  if (length(cell_index) != 1L ||
      is.na(cell_index) ||
      cell_index < 1L ||
      cell_index > length(prepared$cell_id)) {
    .fire_stop('cell_index must identify one prepared active cell.')
  }

  if (!is.numeric(current_time) ||
      length(current_time) != 1L ||
      is.na(current_time) ||
      !is.finite(current_time)) {
    .fire_stop('current_time must be one finite numeric value.')
  }

  if (!is.numeric(timestep) ||
      length(timestep) != 1L ||
      is.na(timestep) ||
      !is.finite(timestep) ||
      timestep <= 0) {
    .fire_stop('timestep must be one finite value > 0.')
  }

  region_index <- prepared$region_index[cell_index]

  # FireEvent.Initiate:
  # annual IgnitionProbability is multiplied directly by timestep.
  # It is intentionally NOT clamped to 1.
  ignition_probability <-
    prepared$ignition_probability[cell_index] * timestep

  u_ignition <- .basefire_uniform_one(uniform)

  result <- list(
    initiated = FALSE,
    cell_index = cell_index,
    cell_id = prepared$cell_id[cell_index],
    region_index = region_index,
    ignition_probability = ignition_probability,
    u_ignition = u_ignition,
    time_since_last_fire = NA_real_,
    init_spread_probability = NA_real_,
    u_init_spread = NA_real_,
    severity = NA_real_,
    draws = 1L,
    failed_gate = 'ignition'
  )

  # Preserve C# short-circuit evaluation and RNG consumption.
  if (!(u_ignition <= ignition_probability)) {
    class(result) <- 'landisR_fire_basefire_ignition_result'
    return(result)
  }

  time_since_last_fire <-
    current_time - prepared$last_fire_year[cell_index]

  init_spread_probability <- .basefire_compute_init_spread_prob(
    time_since_last_fire = time_since_last_fire,
    fire_spread_age = prepared$fire_spread_age[cell_index]
  )

  u_init_spread <- .basefire_uniform_one(uniform)

  result$time_since_last_fire <- time_since_last_fire
  result$init_spread_probability <- init_spread_probability
  result$u_init_spread <- u_init_spread
  result$draws <- 2L
  result$failed_gate <- 'init_spread'

  # A comparison with NaN is FALSE, matching C# double NaN behavior.
  if (!isTRUE(u_init_spread <= init_spread_probability)) {
    class(result) <- 'landisR_fire_basefire_ignition_result'
    return(result)
  }

  severity <- severity_fun(
    cell_index = cell_index,
    current_time = current_time
  )

  if (!is.numeric(severity) ||
      length(severity) != 1L ||
      is.na(severity) ||
      !is.finite(severity)) {
    .fire_stop('severity_fun must return one finite numeric severity.')
  }

  result$severity <- as.numeric(severity)
  result$failed_gate <- 'severity'

  if (!(severity > 0)) {
    class(result) <- 'landisR_fire_basefire_ignition_result'
    return(result)
  }

  result$initiated <- TRUE
  result$failed_gate <- NA_character_

  class(result) <- 'landisR_fire_basefire_ignition_result'

  if (isTRUE(diagnostics)) {
    return(result)
  }

  result
}

.basefire_ignition <- .wrap_component(
  fun = function(
      cell_index,
      current_time,
      timestep,
      severity_fun,
      prepared,
      uniform = stats::runif,
      diagnostics = FALSE,
      ...) {

    .basefire_try_initiate(
      cell_index = cell_index,
      current_time = current_time,
      timestep = timestep,
      prepared = prepared,
      severity_fun = severity_fun,
      uniform = uniform,
      diagnostics = diagnostics
    )
  },
  class_name = 'ignition',
  prepare = .basefire_prepare_ignition,
  metadata = list(
    profile = 'basefire',
    reference = 'FireEvent.Initiate'
  )
)
