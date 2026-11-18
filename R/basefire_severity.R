.basefire_align_optional_site_vector <- function(
    value,
    cell_id,
    field_name) {

  if (is.null(value)) {
    return(NULL)
  }

  n <- length(cell_id)

  if (is.data.frame(value)) {
    required <- c('cell_id', field_name)
    missing <- setdiff(required, names(value))

    if (length(missing)) {
      .fire_stop(
        field_name,
        ' data.frame is missing columns: ',
        paste(missing, collapse = ', '),
        '.'
      )
    }

    idx <- match(cell_id, value$cell_id)

    if (anyNA(idx)) {
      .fire_stop(
        field_name,
        ' data.frame does not contain every active cell.'
      )
    }

    value <- value[[field_name]][idx]
  } else if (!is.null(names(value)) && all(nzchar(names(value)))) {
    idx <- match(as.character(cell_id), names(value))

    if (anyNA(idx)) {
      .fire_stop(
        field_name,
        ' named vector does not contain every active cell.'
      )
    }

    value <- value[idx]
  }

  if (!is.numeric(value) || length(value) != n) {
    .fire_stop(
      field_name,
      ' must be numeric and aligned to all active cells.'
    )
  }

  if (anyNA(value) || any(!is.finite(value))) {
    .fire_stop(
      field_name,
      ' cannot contain NA or non-finite values.'
    )
  }

  as.numeric(value)
}

.basefire_prepare_severity <- function(
    landscape,
    fire_state,
    topology,
    context,
    model) {

  cache <- context$basefire_cache

  if (is.null(cache) ||
      !inherits(cache, 'landisR_fire_basefire_parameter_cache')) {
    .fire_stop(
      'Base Fire severity preparation requires context$basefire_cache.'
    )
  }

  if (is.null(topology)) {
    topology <- fire_grid_topology(landscape)
  }

  .validate_fire_preparation_topology(topology)

  region_index <- .basefire_resolve_region_index(
    landscape = landscape,
    topology = topology,
    cache = cache
  )

  ignition_cache <- .basefire_prepare_ignition(
    landscape = landscape,
    fire_state = fire_state,
    topology = topology,
    context = context,
    model = model
  )

  last_wind_year <- NULL

  if (!is.null(fire_state) &&
      'last_wind_year' %in% names(fire_state$cells)) {

    state_idx <- match(
      topology$cell_id,
      fire_state$cells$cell_id
    )

    if (anyNA(state_idx)) {
      .fire_stop(
        'fire_state is missing at least one active topology cell.'
      )
    }

    last_wind_year <- .basefire_align_optional_site_vector(
      fire_state$cells$last_wind_year[state_idx],
      topology$cell_id,
      'last_wind_year'
    )
  } else if (!is.null(context$time_of_last_wind)) {
    last_wind_year <- .basefire_align_optional_site_vector(
      context$time_of_last_wind,
      topology$cell_id,
      'last_wind_year'
    )
  }

  x <- list(
    cell_id = topology$cell_id,
    region_index = region_index,
    last_fire_year = ignition_cache$last_fire_year,
    has_wind_state = !is.null(last_wind_year),
    last_wind_year = last_wind_year,
    fuel_curve = cache$fuel_curve,
    wind_curve = cache$wind_curve
  )

  class(x) <- 'landisR_fire_basefire_severity_cache'
  .basefire_validate_severity_cache(x)
  x
}

.basefire_validate_severity_cache <- function(x) {
  if (!inherits(x, 'landisR_fire_basefire_severity_cache')) {
    .fire_stop(
      'prepared must inherit from landisR_fire_basefire_severity_cache.'
    )
  }

  n <- length(x$cell_id)

  if (length(x$region_index) != n ||
      length(x$last_fire_year) != n) {
    .fire_stop('Invalid Base Fire severity cache lengths.')
  }

  if (!is.logical(x$has_wind_state) ||
      length(x$has_wind_state) != 1L ||
      is.na(x$has_wind_state)) {
    .fire_stop('Invalid Base Fire has_wind_state flag.')
  }

  if (x$has_wind_state) {
    if (!is.numeric(x$last_wind_year) ||
        length(x$last_wind_year) != n ||
        anyNA(x$last_wind_year) ||
        any(!is.finite(x$last_wind_year))) {
      .fire_stop('Invalid Base Fire last_wind_year vector.')
    }
  } else if (!is.null(x$last_wind_year)) {
    .fire_stop(
      'last_wind_year must be NULL when has_wind_state is FALSE.'
    )
  }

  if (!is.matrix(x$fuel_curve) ||
      ncol(x$fuel_curve) != 5L ||
      !is.matrix(x$wind_curve) ||
      ncol(x$wind_curve) != 5L) {
    .fire_stop('Base Fire fuel/wind curves must be region-by-5 matrices.')
  }

  invisible(TRUE)
}

.basefire_curve_severity_fuel <- function(time_since_last_fire, thresholds) {
  if (!is.numeric(time_since_last_fire) ||
      length(time_since_last_fire) != 1L ||
      is.na(time_since_last_fire) ||
      !is.finite(time_since_last_fire)) {
    .fire_stop(
      'time_since_last_fire must be one finite numeric value.'
    )
  }

  if (!is.numeric(thresholds) ||
      length(thresholds) != 5L ||
      anyNA(thresholds)) {
    .fire_stop('Fuel thresholds must contain five numeric values.')
  }

  severity <- 0L

  for (j in seq_len(5L)) {
    if (thresholds[j] != -1 &&
        time_since_last_fire >= thresholds[j]) {
      severity <- j
    }
  }

  as.integer(severity)
}

.basefire_curve_severity_wind <- function(time_since_last_wind, thresholds) {
  if (!is.numeric(time_since_last_wind) ||
      length(time_since_last_wind) != 1L ||
      is.na(time_since_last_wind) ||
      !is.finite(time_since_last_wind)) {
    .fire_stop(
      'time_since_last_wind must be one finite numeric value.'
    )
  }

  if (!is.numeric(thresholds) ||
      length(thresholds) != 5L ||
      anyNA(thresholds)) {
    .fire_stop('Wind thresholds must contain five numeric values.')
  }

  severity <- 0L

  for (j in seq_len(5L)) {
    if (thresholds[j] != -1 &&
        time_since_last_wind <= thresholds[j]) {
      severity <- j
    }
  }

  as.integer(severity)
}

.basefire_calc_severity <- function(
    cell_index,
    current_time,
    prepared,
    last_fire_year = NULL,
    last_wind_year = NULL,
    diagnostics = FALSE) {

  .basefire_validate_severity_cache(prepared)

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

  if (is.null(last_fire_year)) {
    last_fire_year <- prepared$last_fire_year
  }

  if (!is.numeric(last_fire_year) ||
      length(last_fire_year) != length(prepared$cell_id)) {
    .fire_stop(
      'last_fire_year must match the prepared active-cell length.'
    )
  }

  region <- prepared$region_index[cell_index]
  time_since_last_fire <- current_time - last_fire_year[cell_index]

  fuel_severity <- .basefire_curve_severity_fuel(
    time_since_last_fire = time_since_last_fire,
    thresholds = prepared$fuel_curve[region, ]
  )

  wind_severity <- 0L
  time_since_last_wind <- NA_real_

  has_wind <- prepared$has_wind_state || !is.null(last_wind_year)

  if (has_wind) {
    if (is.null(last_wind_year)) {
      last_wind_year <- prepared$last_wind_year
    }

    if (!is.numeric(last_wind_year) ||
        length(last_wind_year) != length(prepared$cell_id)) {
      .fire_stop(
        'last_wind_year must match the prepared active-cell length.'
      )
    }

    time_since_last_wind <- current_time - last_wind_year[cell_index]

    wind_severity <- .basefire_curve_severity_wind(
      time_since_last_wind = time_since_last_wind,
      thresholds = prepared$wind_curve[region, ]
    )
  }

  severity <- max(fuel_severity, wind_severity)

  if (!isTRUE(diagnostics)) {
    return(as.integer(severity))
  }

  out <- list(
    cell_index = cell_index,
    cell_id = prepared$cell_id[cell_index],
    region_index = region,
    time_since_last_fire = time_since_last_fire,
    fuel_severity = fuel_severity,
    has_wind_state = has_wind,
    time_since_last_wind = time_since_last_wind,
    wind_severity = wind_severity,
    severity = as.integer(severity)
  )

  class(out) <- 'landisR_fire_basefire_severity_result'
  out
}

.basefire_severity <- .wrap_component(
  fun = function(
      cell_index,
      current_time,
      prepared,
      last_fire_year = NULL,
      last_wind_year = NULL,
      diagnostics = FALSE,
      ...) {

    .basefire_calc_severity(
      cell_index = cell_index,
      current_time = current_time,
      prepared = prepared,
      last_fire_year = last_fire_year,
      last_wind_year = last_wind_year,
      diagnostics = diagnostics
    )
  },
  class_name = 'severity',
  prepare = .basefire_prepare_severity,
  metadata = list(
    profile = 'basefire',
    reference = 'FireEvent.CalcSeverity'
  )
)
