.basefire_float32 <- function(x) {
  if (!is.numeric(x)) {
    .fire_stop('x must be numeric for float32 conversion.')
  }

  if (!length(x)) {
    return(numeric())
  }

  raw <- writeBin(
    as.numeric(x),
    raw(),
    size = 4L,
    endian = .Platform$endian
  )

  readBin(
    raw,
    what = 'numeric',
    n = length(x),
    size = 4L,
    endian = .Platform$endian
  )
}

.basefire_normalize_species_longevity <- function(
    species_longevity = NULL,
    species = NULL) {

  if (!is.null(species_longevity)) {
    if (!is.numeric(species_longevity) ||
        is.null(names(species_longevity)) ||
        any(!nzchar(names(species_longevity)))) {
      .fire_stop(
        'species_longevity must be a named numeric vector.'
      )
    }

    if (anyDuplicated(names(species_longevity))) {
      .fire_stop('species_longevity names must be unique.')
    }

    longevity <- as.numeric(species_longevity)
    names(longevity) <- names(species_longevity)
  } else {
    if (is.null(species) || !is.data.frame(species)) {
      .fire_stop(
        'Base Fire mortality preparation requires either ',
        'context$species_longevity or context$species.'
      )
    }

    .basefire_require_columns(
      species,
      c('species_code', 'longevity'),
      'context$species'
    )

    species_code <- .basefire_nonempty_character(
      species$species_code,
      'context$species$species_code'
    )

    if (anyDuplicated(species_code)) {
      .fire_stop('context$species species_code values must be unique.')
    }

    longevity <- as.numeric(species$longevity)

    if (anyNA(longevity) || any(!is.finite(longevity))) {
      .fire_stop(
        'context$species$longevity must contain finite numeric values.'
      )
    }

    names(longevity) <- species_code
  }

  if (anyNA(longevity) ||
      any(!is.finite(longevity)) ||
      any(longevity <= 0)) {
    .fire_stop('Species longevity values must be finite and > 0.')
  }

  longevity
}

.basefire_prepare_mortality <- function(
    landscape,
    fire_state,
    topology,
    context,
    model) {

  cache <- context$basefire_cache

  if (is.null(cache) ||
      !inherits(cache, 'landisR_fire_basefire_parameter_cache')) {
    .fire_stop(
      'Base Fire mortality preparation requires context$basefire_cache.'
    )
  }

  longevity <- .basefire_normalize_species_longevity(
    species_longevity = context$species_longevity,
    species = context$species
  )

  required_species <- names(cache$species_fire_tolerance)
  missing_longevity <- setdiff(required_species, names(longevity))

  if (length(missing_longevity)) {
    .fire_stop(
      'Missing longevity for Base Fire species: ',
      paste(missing_longevity, collapse = ', '),
      '.'
    )
  }

  # Preserve the Base Fire species-file order in compact aligned vectors.
  longevity <- longevity[required_species]

  x <- list(
    species_code = required_species,
    species_fire_tolerance = as.integer(
      unname(cache$species_fire_tolerance[required_species])
    ),
    species_longevity = as.numeric(unname(longevity)),
    species_index = setNames(
      seq_along(required_species),
      required_species
    ),
    damage_max_age = as.numeric(cache$damage_max_age),
    damage_severity_tolerance_diff = as.integer(
      cache$damage_severity_tolerance_diff
    )
  )

  class(x) <- 'landisR_fire_basefire_mortality_cache'
  .basefire_validate_mortality_cache(x)
  x
}

.basefire_validate_mortality_cache <- function(x) {
  if (!inherits(x, 'landisR_fire_basefire_mortality_cache')) {
    .fire_stop(
      'prepared must inherit from landisR_fire_basefire_mortality_cache.'
    )
  }

  n_species <- length(x$species_code)

  if (!n_species ||
      length(x$species_fire_tolerance) != n_species ||
      length(x$species_longevity) != n_species) {
    .fire_stop('Invalid Base Fire mortality species cache.')
  }

  if (anyDuplicated(x$species_code) ||
      anyNA(x$species_fire_tolerance) ||
      anyNA(x$species_longevity) ||
      any(x$species_longevity <= 0)) {
    .fire_stop('Invalid Base Fire mortality species values.')
  }

  n_damage <- length(x$damage_max_age)

  if (!n_damage ||
      length(x$damage_severity_tolerance_diff) != n_damage) {
    .fire_stop('Invalid Base Fire mortality damage cache.')
  }

  if (anyNA(x$damage_max_age) ||
      anyNA(x$damage_severity_tolerance_diff)) {
    .fire_stop('Base Fire mortality damage values cannot contain NA.')
  }

  invisible(TRUE)
}

.basefire_mortality_species_index <- function(
    species_code,
    prepared) {

  if (!is.character(species_code) ||
      length(species_code) != 1L ||
      is.na(species_code) ||
      !nzchar(species_code)) {
    .fire_stop('species_code must be one non-empty character value.')
  }

  idx <- unname(prepared$species_index[species_code])

  if (!length(idx) || is.na(idx)) {
    .fire_stop(
      'Unknown Base Fire species code: ',
      species_code,
      '.'
    )
  }

  as.integer(idx)
}

.basefire_cohort_mortality <- function(
    severity,
    species_code,
    cohort_age,
    prepared,
    diagnostics = FALSE) {

  .basefire_validate_mortality_cache(prepared)

  if (!is.numeric(severity) ||
      length(severity) != 1L ||
      is.na(severity) ||
      !is.finite(severity) ||
      severity < 0 ||
      abs(severity - round(severity)) > sqrt(.Machine$double.eps)) {
    .fire_stop('severity must be one non-negative integer-like value.')
  }

  severity <- as.integer(round(severity))

  if (!is.numeric(cohort_age) ||
      length(cohort_age) != 1L ||
      is.na(cohort_age) ||
      !is.finite(cohort_age) ||
      cohort_age < 0) {
    .fire_stop('cohort_age must be one finite numeric value >= 0.')
  }

  # Exact first branch of ReduceOrKillMarkedCohort():
  # severity 5 kills every cohort without consulting tolerance,
  # longevity or the damage table.
  if (severity == 5L) {
    out <- list(
      killed = TRUE,
      removal_fraction = 1.0,
      severity = severity,
      species_code = species_code,
      cohort_age = cohort_age,
      fire_tolerance = NA_integer_,
      severity_tolerance_difference = NA_integer_,
      longevity = NA_real_,
      age_as_percent = NA_real_,
      matched_damage_index = NA_integer_,
      matched_damage_max_age = NA_real_,
      matched_damage_severity_tolerance_diff = NA_integer_,
      reason = 'severity5'
    )

    class(out) <- 'landisR_fire_basefire_mortality_result'

    if (isTRUE(diagnostics)) out else 1.0
  } else {
    species_index <- .basefire_mortality_species_index(
      species_code = species_code,
      prepared = prepared
    )

    fire_tolerance <- prepared$species_fire_tolerance[species_index]
    longevity <- prepared$species_longevity[species_index]

    # C#:
    # float ageAsPercent =
    #   (float) cohort.Data.Age / (float) cohort.Species.Longevity;
    #
    # Reproduce the single-precision casts and result.
    age_float <- .basefire_float32(cohort_age)
    longevity_float <- .basefire_float32(longevity)
    age_as_percent <- .basefire_float32(
      age_float / longevity_float
    )

    severity_tolerance_difference <- severity - fire_tolerance

    kill_cohort <- FALSE
    matched <- NA_integer_

    # The current parser appends FireDamages in increasing input order and
    # PlugIn passes that list directly to FireEvent.Initialize(). Although
    # an old source comment says "read table backwards", the executed foreach
    # iterates this list forward. Preserve executed behavior.
    for (i in seq_along(prepared$damage_max_age)) {
      if (
        severity_tolerance_difference >=
          prepared$damage_severity_tolerance_diff[i]
      ) {
        max_age_float <- .basefire_float32(
          prepared$damage_max_age[i]
        )

        if (max_age_float >= age_as_percent) {
          kill_cohort <- TRUE
          matched <- i
          break
        }
      }
    }

    out <- list(
      killed = kill_cohort,
      removal_fraction = if (kill_cohort) 1.0 else 0.0,
      severity = severity,
      species_code = species_code,
      cohort_age = cohort_age,
      fire_tolerance = fire_tolerance,
      severity_tolerance_difference =
        severity_tolerance_difference,
      longevity = longevity,
      age_as_percent = age_as_percent,
      matched_damage_index = matched,
      matched_damage_max_age =
        if (is.na(matched)) NA_real_
        else prepared$damage_max_age[matched],
      matched_damage_severity_tolerance_diff =
        if (is.na(matched)) NA_integer_
        else prepared$damage_severity_tolerance_diff[matched],
      reason = if (kill_cohort) 'damage_table' else 'survives'
    )

    class(out) <- 'landisR_fire_basefire_mortality_result'

    if (isTRUE(diagnostics)) {
      out
    } else if (kill_cohort) {
      1.0
    } else {
      0.0
    }
  }
}

.basefire_apply_mortality <- function(
    cohorts,
    severity,
    prepared,
    keep_diagnostics = FALSE) {

  .basefire_validate_mortality_cache(prepared)

  if (!is.data.frame(cohorts)) {
    .fire_stop('cohorts must be a data.frame.')
  }

  .basefire_require_columns(
    cohorts,
    c('species_code', 'age'),
    'cohorts'
  )

  n <- nrow(cohorts)

  if (!n) {
    out <- list(
      cohorts_killed = 0L,
      kill = logical(),
      removal_fraction = numeric(),
      survivors = cohorts,
      killed = cohorts,
      diagnostics = if (isTRUE(keep_diagnostics)) list() else NULL
    )

    class(out) <- 'landisR_fire_basefire_site_mortality_result'
    return(out)
  }

  if (!is.numeric(severity) ||
      length(severity) != 1L ||
      is.na(severity) ||
      !is.finite(severity)) {
    .fire_stop('severity must be one finite numeric value.')
  }

  kill <- logical(n)
  removal_fraction <- numeric(n)
  diagnostic_list <- if (isTRUE(keep_diagnostics)) vector('list', n) else NULL

  for (i in seq_len(n)) {
    result <- .basefire_cohort_mortality(
      severity = severity,
      species_code = as.character(cohorts$species_code[i]),
      cohort_age = cohorts$age[i],
      prepared = prepared,
      diagnostics = TRUE
    )

    kill[i] <- result$killed
    removal_fraction[i] <- result$removal_fraction

    if (isTRUE(keep_diagnostics)) {
      diagnostic_list[[i]] <- result
    }
  }

  out <- list(
    cohorts_killed = as.integer(sum(kill)),
    kill = kill,
    removal_fraction = removal_fraction,
    survivors = cohorts[!kill, , drop = FALSE],
    killed = cohorts[kill, , drop = FALSE],
    diagnostics = diagnostic_list
  )

  rownames(out$survivors) <- NULL
  rownames(out$killed) <- NULL

  class(out) <- 'landisR_fire_basefire_site_mortality_result'
  out
}

.basefire_mortality <- .wrap_component(
  fun = function(
      severity,
      species_code,
      cohort_age,
      prepared,
      diagnostics = FALSE,
      ...) {

    .basefire_cohort_mortality(
      severity = severity,
      species_code = species_code,
      cohort_age = cohort_age,
      prepared = prepared,
      diagnostics = diagnostics
    )
  },
  class_name = 'mortality',
  prepare = .basefire_prepare_mortality,
  metadata = list(
    profile = 'basefire',
    reference = 'FireEvent.ReduceOrKillMarkedCohort'
  )
)
