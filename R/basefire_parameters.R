.basefire_empty_dynamic_regions <- function() {
  data.frame(
    year = integer(),
    map_name = character(),
    stringsAsFactors = FALSE
  )
}

.basefire_empty_curves <- function() {
  data.frame(
    region_name = character(),
    severity1 = integer(),
    severity2 = integer(),
    severity3 = integer(),
    severity4 = integer(),
    severity5 = integer(),
    stringsAsFactors = FALSE
  )
}

.basefire_require_columns <- function(x, required, object_name) {
  if (!is.data.frame(x)) {
    .fire_stop(object_name, ' must be a data.frame.')
  }

  missing <- setdiff(required, names(x))
  if (length(missing)) {
    .fire_stop(
      object_name,
      ' is missing required columns: ',
      paste(missing, collapse = ', '),
      '.'
    )
  }

  invisible(TRUE)
}

.basefire_integer_like <- function(x, name) {
  if (!is.numeric(x) || any(!is.finite(x))) {
    .fire_stop(name, ' must contain finite numeric values.')
  }

  rounded <- round(x)
  tol <- sqrt(.Machine$double.eps)

  if (any(abs(x - rounded) > tol)) {
    .fire_stop(name, ' must contain integer-like values.')
  }

  as.integer(rounded)
}

.basefire_nonempty_character <- function(x, name) {
  if (!is.character(x) || anyNA(x) || any(!nzchar(trimws(x)))) {
    .fire_stop(name, ' must contain non-empty character values.')
  }

  trimws(x)
}

.basefire_normalize_regions <- function(regions) {
  required <- c(
    'name',
    'map_code',
    'mean_size',
    'min_size',
    'max_size',
    'ignition_probability',
    'fire_spread_age'
  )

  .basefire_require_columns(regions, required, 'regions')

  name <- .basefire_nonempty_character(regions$name, 'regions$name')
  map_code <- .basefire_integer_like(regions$map_code, 'regions$map_code')
  fire_spread_age <- .basefire_integer_like(
    regions$fire_spread_age,
    'regions$fire_spread_age'
  )

  numeric_cols <- c(
    'mean_size',
    'min_size',
    'max_size',
    'ignition_probability'
  )

  for (nm in numeric_cols) {
    if (!is.numeric(regions[[nm]]) || any(!is.finite(regions[[nm]]))) {
      .fire_stop('regions$', nm, ' must contain finite numeric values.')
    }
  }

  if (!length(name)) {
    .fire_stop('At least one Base Fire region is required.')
  }
  if (anyDuplicated(name)) {
    .fire_stop('Base Fire region names must be unique.')
  }
  if (any(map_code < 0L | map_code > 65535L)) {
    .fire_stop('Base Fire map_code values must be between 0 and 65535.')
  }
  if (anyDuplicated(map_code)) {
    .fire_stop('Base Fire map_code values must be unique.')
  }
  if (any(regions$mean_size < 0) ||
      any(regions$min_size < 0) ||
      any(regions$max_size < 0)) {
    .fire_stop('Base Fire fire-size parameters must be >= 0.')
  }
  if (any(regions$ignition_probability < 0 |
          regions$ignition_probability > 1)) {
    .fire_stop('Base Fire ignition_probability must be between 0 and 1.')
  }
  if (any(fire_spread_age < 0L)) {
    .fire_stop('Base Fire fire_spread_age must be >= 0.')
  }

  # LANDIS-II Base Fire does not enforce an ordering relationship among
  # MinSize, MeanSize and MaxSize in FireRegion/InputParameterParser.
  data.frame(
    region_index = seq_along(name),
    basefire_index = seq_along(name) - 1L,
    name = name,
    map_code = map_code,
    mean_size = as.numeric(regions$mean_size),
    min_size = as.numeric(regions$min_size),
    max_size = as.numeric(regions$max_size),
    ignition_probability = as.numeric(regions$ignition_probability),
    fire_spread_age = fire_spread_age,
    stringsAsFactors = FALSE
  )
}

.basefire_normalize_dynamic_regions <- function(dynamic_regions) {
  if (is.null(dynamic_regions)) {
    dynamic_regions <- .basefire_empty_dynamic_regions()
  }

  .basefire_require_columns(
    dynamic_regions,
    c('year', 'map_name'),
    'dynamic_regions'
  )

  if (!nrow(dynamic_regions)) {
    return(.basefire_empty_dynamic_regions())
  }

  year <- .basefire_integer_like(dynamic_regions$year, 'dynamic_regions$year')
  map_name <- .basefire_nonempty_character(
    dynamic_regions$map_name,
    'dynamic_regions$map_name'
  )

  if (any(year <= 0L)) {
    .fire_stop('Dynamic Base Fire region years must be > 0.')
  }
  if (any(diff(year) <= 0L)) {
    .fire_stop('Dynamic Base Fire region years must be strictly increasing.')
  }

  data.frame(
    year = year,
    map_name = map_name,
    stringsAsFactors = FALSE
  )
}

.basefire_normalize_curve <- function(curves, regions, type = c('fuel', 'wind')) {
  type <- match.arg(type)

  if (is.null(curves)) {
    curves <- .basefire_empty_curves()
  }

  required <- c(
    'region_name',
    'severity1',
    'severity2',
    'severity3',
    'severity4',
    'severity5'
  )

  .basefire_require_columns(curves, required, paste0(type, '_curves'))

  if (!nrow(curves)) {
    return(.basefire_empty_curves())
  }

  region_name <- .basefire_nonempty_character(
    curves$region_name,
    paste0(type, '_curves$region_name')
  )

  if (anyDuplicated(region_name)) {
    .fire_stop(type, ' curve region names must be unique.')
  }

  unknown <- setdiff(region_name, regions$name)
  if (length(unknown)) {
    .fire_stop(
      type,
      ' curve references unknown Base Fire region(s): ',
      paste(unknown, collapse = ', '),
      '.'
    )
  }

  severity <- matrix(
    0L,
    nrow = nrow(curves),
    ncol = 5L,
    dimnames = list(NULL, paste0('severity', 1:5))
  )

  for (j in seq_len(5L)) {
    severity[, j] <- .basefire_integer_like(
      curves[[paste0('severity', j)]],
      paste0(type, '_curves$severity', j)
    )
  }

  if (any(severity < -1L)) {
    .fire_stop(type, ' curve severity ages must be >= -1.')
  }

  if (nrow(severity)) {
    for (i in seq_len(nrow(severity))) {
      v <- severity[i, ]

      for (j in 2:5) {
        if (type == 'fuel') {
          # Exact LANDIS-II FuelCurve setter behavior:
          # a non-negative threshold must be strictly greater than the
          # preceding threshold.
          if (v[j] > -1L && v[j] <= v[j - 1L]) {
            .fire_stop(
              'Invalid Base Fire fuel curve for region ',
              region_name[i],
              ': severity',
              j,
              ' must be -1 or > severity',
              j - 1L,
              '.'
            )
          }
        } else {
          # Exact LANDIS-II WindCurve setter behavior:
          # equality is accepted; only a decrease is rejected.
          if (v[j] != -1L && v[j] < v[j - 1L]) {
            .fire_stop(
              'Invalid Base Fire wind curve for region ',
              region_name[i],
              ': severity',
              j,
              ' must be -1 or >= severity',
              j - 1L,
              '.'
            )
          }
        }
      }
    }
  }

  data.frame(
    region_name = region_name,
    severity1 = severity[, 1L],
    severity2 = severity[, 2L],
    severity3 = severity[, 3L],
    severity4 = severity[, 4L],
    severity5 = severity[, 5L],
    stringsAsFactors = FALSE
  )
}

.basefire_normalize_damage <- function(damage) {
  .basefire_require_columns(
    damage,
    c('max_age', 'severity_tolerance_diff'),
    'damage'
  )

  if (!nrow(damage)) {
    .fire_stop('At least one Base Fire damage class is required.')
  }

  if (!is.numeric(damage$max_age) || any(!is.finite(damage$max_age))) {
    .fire_stop('damage$max_age must contain finite numeric values.')
  }

  max_age <- as.numeric(damage$max_age)
  severity_tolerance_diff <- .basefire_integer_like(
    damage$severity_tolerance_diff,
    'damage$severity_tolerance_diff'
  )

  if (any(max_age <= 0 | max_age > 1)) {
    .fire_stop('Base Fire damage max_age must be > 0 and <= 1.')
  }
  if (any(diff(max_age) <= 0)) {
    .fire_stop('Base Fire damage max_age must be strictly increasing.')
  }

  # DamageTable itself allows -4..4, while InputParameterParser starts the
  # ordered table at previousNumber=-4 and requires every next value to be
  # strictly greater. Therefore an input-table first value must be > -4.
  if (any(severity_tolerance_diff <= -4L |
          severity_tolerance_diff > 4L)) {
    .fire_stop(
      'Base Fire damage severity_tolerance_diff must be > -4 and <= 4.'
    )
  }
  if (any(diff(severity_tolerance_diff) <= 0L)) {
    .fire_stop(
      'Base Fire damage severity_tolerance_diff must be strictly increasing.'
    )
  }

  data.frame(
    max_age = max_age,
    severity_tolerance_diff = severity_tolerance_diff,
    stringsAsFactors = FALSE
  )
}

.basefire_normalize_species <- function(species) {
  .basefire_require_columns(
    species,
    c('species_code', 'fire_tolerance'),
    'species'
  )

  species_code <- .basefire_nonempty_character(
    species$species_code,
    'species$species_code'
  )

  if (anyDuplicated(species_code)) {
    .fire_stop('Base Fire species codes must be unique.')
  }

  fire_tolerance <- .basefire_integer_like(
    species$fire_tolerance,
    'species$fire_tolerance'
  )

  data.frame(
    species_code = species_code,
    fire_tolerance = fire_tolerance,
    stringsAsFactors = FALSE
  )
}

.basefire_parameters <- function(
    timestep,
    species,
    regions,
    initial_fire_regions,
    dynamic_regions = NULL,
    fuel_curves = NULL,
    wind_curves = NULL,
    damage,
    map_names_template,
    log_file,
    summary_log_file,
    source_file = NULL,
    species_csv_file = NULL) {

  timestep <- .basefire_integer_like(timestep, 'timestep')
  if (length(timestep) != 1L || is.na(timestep) || timestep <= 0L) {
    .fire_stop('Base Fire timestep must be one integer > 0.')
  }

  .assert_scalar_character(initial_fire_regions, 'initial_fire_regions')
  .assert_scalar_character(map_names_template, 'map_names_template')
  .assert_scalar_character(log_file, 'log_file')
  .assert_scalar_character(summary_log_file, 'summary_log_file')

  regions <- .basefire_normalize_regions(regions)
  dynamic_regions <- .basefire_normalize_dynamic_regions(dynamic_regions)
  species <- .basefire_normalize_species(species)
  fuel_curves <- .basefire_normalize_curve(
    fuel_curves,
    regions,
    type = 'fuel'
  )
  wind_curves <- .basefire_normalize_curve(
    wind_curves,
    regions,
    type = 'wind'
  )
  damage <- .basefire_normalize_damage(damage)

  x <- list(
    timestep = timestep,
    species = species,
    regions = regions,
    initial_fire_regions = initial_fire_regions,
    dynamic_regions = dynamic_regions,
    fuel_curves = fuel_curves,
    wind_curves = wind_curves,
    damage = damage,
    map_names_template = map_names_template,
    log_file = log_file,
    summary_log_file = summary_log_file,
    source_file = source_file,
    species_csv_file = species_csv_file
  )

  class(x) <- 'landisR_fire_basefire_parameters'
  .basefire_validate_parameters(x)
  x
}

.basefire_validate_parameters <- function(x) {
  if (!inherits(x, 'landisR_fire_basefire_parameters')) {
    .fire_stop('x must inherit from landisR_fire_basefire_parameters.')
  }

  .basefire_normalize_regions(x$regions)
  .basefire_normalize_dynamic_regions(x$dynamic_regions)
  .basefire_normalize_species(x$species)
  .basefire_normalize_curve(x$fuel_curves, x$regions, type = 'fuel')
  .basefire_normalize_curve(x$wind_curves, x$regions, type = 'wind')
  .basefire_normalize_damage(x$damage)

  invisible(TRUE)
}

.basefire_prepare_parameter_cache <- function(parameters) {
  .basefire_validate_parameters(parameters)

  regions <- parameters$regions
  n_regions <- nrow(regions)

  map_code_to_region <- integer(65536L)
  map_code_to_region[regions$map_code + 1L] <- regions$region_index

  region_by_name <- setNames(regions$region_index, regions$name)

  curve_names <- paste0('severity', 1:5)

  fuel_curve <- matrix(
    0L,
    nrow = n_regions,
    ncol = 5L,
    dimnames = list(regions$name, curve_names)
  )
  wind_curve <- fuel_curve

  if (nrow(parameters$fuel_curves)) {
    idx <- match(parameters$fuel_curves$region_name, regions$name)
    fuel_curve[idx, ] <- as.matrix(
      parameters$fuel_curves[, curve_names, drop = FALSE]
    )
  }

  if (nrow(parameters$wind_curves)) {
    idx <- match(parameters$wind_curves$region_name, regions$name)
    wind_curve[idx, ] <- as.matrix(
      parameters$wind_curves[, curve_names, drop = FALSE]
    )
  }

  storage.mode(fuel_curve) <- 'integer'
  storage.mode(wind_curve) <- 'integer'

  x <- list(
    n_regions = n_regions,
    region_name = regions$name,
    region_index = regions$region_index,
    basefire_index = regions$basefire_index,
    map_code = regions$map_code,
    map_code_to_region = map_code_to_region,
    region_by_name = region_by_name,
    mean_size = regions$mean_size,
    min_size = regions$min_size,
    max_size = regions$max_size,
    ignition_probability = regions$ignition_probability,
    fire_spread_age = regions$fire_spread_age,
    fuel_curve = fuel_curve,
    wind_curve = wind_curve,
    dynamic_year = parameters$dynamic_regions$year,
    dynamic_map_name = parameters$dynamic_regions$map_name,
    damage_max_age = parameters$damage$max_age,
    damage_severity_tolerance_diff =
      parameters$damage$severity_tolerance_diff,
    species_fire_tolerance = setNames(
      parameters$species$fire_tolerance,
      parameters$species$species_code
    )
  )

  class(x) <- 'landisR_fire_basefire_parameter_cache'
  x
}
