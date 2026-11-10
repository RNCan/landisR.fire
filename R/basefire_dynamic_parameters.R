# 272 | Dynamic numeric Base Fire region parameters.
# Uses the validated Base Fire preparation path; does not modify the C++ kernel.
# Scheduled changes apply at the first fire step whose simulation year is >= year.

.basefire_272_fields <- c('mean_size', 'min_size', 'max_size',
                         'ignition_probability', 'fire_spread_age')

#' Construct a schedule of regional Base Fire parameter changes
#'
#' @param changes A data.frame with `year` and exactly one of `region_name`
#'   or `map_code`, plus at least one parameter column. `NA` means no change.
#' @return An immutable schedule for `basefire_update_parameters()`.
#' @export
basefire_parameter_schedule <- function(changes) {
  if (inherits(changes, 'landisR_fire_parameter_schedule')) {
    return(changes)
  }
  if (!is.data.frame(changes) || !nrow(changes)) {
    stop('changes must be a nonempty data.frame.', call. = FALSE)
  }
  keys <- intersect(c('region_name', 'map_code'), names(changes))
  if (!'year' %in% names(changes) || length(keys) != 1L) {
    stop('Provide year and exactly one region key: region_name OR map_code.',
         call. = FALSE)
  }
  present <- intersect(.basefire_272_fields, names(changes))
  if (!length(present)) {
    stop('Supply at least one numeric Base Fire region parameter.',
         call. = FALSE)
  }
  unknown <- setdiff(names(changes), c('year', keys, .basefire_272_fields))
  if (length(unknown)) {
    stop('Unsupported schedule columns: ', paste(unknown, collapse = ', '),
         call. = FALSE)
  }
  year <- changes$year
  if (!is.numeric(year) || anyNA(year) || any(!is.finite(year)) ||
      any(year < 0) || any(year != floor(year)) ||
      any(year > .Machine$integer.max)) {
    stop('year must contain finite nonnegative integers.', call. = FALSE)
  }
  key <- keys[[1L]]
  if (identical(key, 'region_name')) {
    names_region <- as.character(changes$region_name)
    if (anyNA(names_region) || any(!nzchar(trimws(names_region)))) {
      stop('region_name must contain nonempty names.', call. = FALSE)
    }
    names_region <- trimws(names_region)
    codes <- rep(NA_integer_, nrow(changes))
    matchkey <- names_region
  } else {
    codes <- changes$map_code
    if (!is.numeric(codes) || anyNA(codes) || any(!is.finite(codes)) ||
        any(codes != floor(codes)) || any(codes < 0 | codes > 65535)) {
      stop('map_code must contain integers in 0..65535.', call. = FALSE)
    }
    codes <- as.integer(codes)
    names_region <- rep(NA_character_, nrow(changes))
    matchkey <- as.character(codes)
  }
  if (anyDuplicated(paste(year, matchkey, sep = ':'))) {
    stop('At most one update per region and year is allowed.', call. = FALSE)
  }
  out <- data.frame(year = as.integer(year),
                    region_name = names_region,
                    map_code = codes,
                    stringsAsFactors = FALSE)
  for (field in .basefire_272_fields) {
    value <- if (field %in% present) changes[[field]] else
      rep(NA_real_, nrow(changes))
    if (!is.numeric(value) || any(!is.finite(value[!is.na(value)]))) {
      stop(field, ' must be numeric and finite where provided.',
           call. = FALSE)
    }
    if (field == 'fire_spread_age' &&
        any(value[!is.na(value)] != floor(value[!is.na(value)]) |
            value[!is.na(value)] > .Machine$integer.max)) {
      stop('fire_spread_age (K) must be integer-valued.', call. = FALSE)
    }
    if (field == 'ignition_probability' &&
        any(value[!is.na(value)] < 0 | value[!is.na(value)] > 1)) {
      stop('ignition_probability must be in [0,1].', call. = FALSE)
    }
    if (field != 'ignition_probability' &&
        any(value[!is.na(value)] < 0)) {
      stop(field, ' cannot be negative.', call. = FALSE)
    }
    out[[field]] <- as.numeric(value)
  }
  if (any(rowSums(!is.na(out[, .basefire_272_fields, drop = FALSE])) == 0L)) {
    stop('Every schedule row must change at least one parameter.',
         call. = FALSE)
  }
  ord <- order(out$year, matchkey)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, 'region_key') <- key
  class(out) <- c('landisR_fire_parameter_schedule', 'data.frame')
  out
}

#' Update the numeric fire-region parameters of a prepared Base Fire object
#'
#' @param prepared A `basefire_prepare()` result or prior updated profile.
#' @param schedule A schedule from `basefire_parameter_schedule()` or its input.
#' @param current_time Simulation year for which the fire step will run.
#' @return A fully re-prepared profile preserving the fire state and dynamic
#'   map reader, with a sparse update history in `context$parameter_change_log`.
#' @export
basefire_update_parameters <- function(prepared, schedule, current_time) {
  .basefire_validate_prepared_profile(prepared)
  schedule <- basefire_parameter_schedule(schedule)
  if (!is.numeric(current_time) || length(current_time) != 1L ||
      is.na(current_time) || !is.finite(current_time) ||
      current_time != floor(current_time) ||
      current_time < 0 || current_time > .Machine$integer.max) {
    stop('current_time must be one nonnegative integer year.', call. = FALSE)
  }
  current_time <- as.integer(current_time)
  regions <- prepared$parameters$regions
  key <- attr(schedule, 'region_key', exact = TRUE)
  if (identical(key, 'region_name')) {
    idx <- match(schedule$region_name, regions$name)
    unknown <- unique(schedule$region_name[is.na(idx)])
  } else if (identical(key, 'map_code')) {
    idx <- match(schedule$map_code, regions$map_code)
    unknown <- unique(schedule$map_code[is.na(idx)])
  } else {
    stop('Invalid parameter schedule region key.', call. = FALSE)
  }
  if (length(unknown)) {
    stop('Unknown Base Fire region(s): ', paste(unknown, collapse=', '),
         call. = FALSE)
  }
  previous_schedule <- prepared$context$parameter_schedule_272
  previous_time <- prepared$context$parameter_last_time_272
  if (!is.null(previous_schedule) && !identical(previous_schedule, schedule)) {
    stop('Schedule cannot change mid-run; create a new prepared profile.',
         call. = FALSE)
  }
  if (is.null(previous_time)) previous_time <- -1L
  if (current_time < previous_time) {
    stop('current_time cannot go backwards for a prepared fire profile.',
         call. = FALSE)
  }
  hit <- which(schedule$year > previous_time & schedule$year <= current_time)
  log <- prepared$context$parameter_change_log
  if (is.null(log)) {
    log <- data.frame(scheduled_year=integer(), applied_year=integer(),
                      region_name=character(), parameter=character(),
                      old_value=numeric(), new_value=numeric(),
                      stringsAsFactors=FALSE)
  }
  for (j in hit) {
    i <- idx[j]
    for (field in .basefire_272_fields) {
      value <- schedule[[field]][j]
      if (is.na(value)) next
      old <- regions[[field]][i]
      if (isTRUE(all.equal(as.numeric(old), value, tolerance=0))) next
      regions[[field]][i] <- if (field == 'fire_spread_age')
        as.integer(value) else as.numeric(value)
      log[nrow(log) + 1L, ] <- list(schedule$year[j], current_time,
                                   regions$name[i], field,
                                   as.numeric(old), as.numeric(value))
    }
  }
  context <- prepared$context
  context$parameter_schedule_272 <- schedule
  context$parameter_last_time_272 <- current_time
  context$parameter_change_log <- log
  # No values changed, so reuse all expensive prepared caches and workspace.
  if (identical(regions, prepared$parameters$regions)) {
    prepared$context <- context
    return(prepared)
  }
  params <- prepared$parameters
  params$regions <- regions
  # Revalidate the *whole* parameter object, including all region constraints.
  .basefire_validate_parameters(params)
  # Build every regional and per-cell cache from the changed parameters.
  # The currently held fire_state is supplied unchanged. This also preserves
  # the dynamic fire-region map reader and time-of-last-wind context.
  profile <- .basefire_profile(params)
  updated <- .basefire_prepare_profile(
    profile=profile,
    landscape=prepared$landscape,
    fire_state=prepared$fire_state,
    context=context,
    topology=NULL
  )
  if (!identical(updated$fire_state, prepared$fire_state)) {
    stop('Unexpected fire_state change during dynamic re-preparation.',
         call. = FALSE)
  }
  updated
}

#' Simulate one Base Fire step after applying scheduled parameter updates
#'
#' @param simulation Current landisR simulation (after succession step).
#' @param prepared Current prepared Base Fire object.
#' @param parameter_schedule Schedule of changes in numeric region parameters.
#' @param ... Passed unchanged to `basefire_step()` including optional output
#'   writing arguments (version 0.1.0.9013 or newer).
#' @return The standard `basefire_step()` result with an updated
#'   `prepared_profile` and unchanged original return field names.
#' @export
basefire_step_dynamic <- function(simulation, prepared,
                                  parameter_schedule, ...) {
  if (is.null(simulation$year)) {
    stop('simulation$year is required for dynamic Base Fire updates.',
         call. = FALSE)
  }
  prepared <- basefire_update_parameters(
    prepared=prepared, schedule=parameter_schedule,
    current_time=simulation$year
  )
  basefire_step(simulation=simulation, prepared=prepared, ...)
}
