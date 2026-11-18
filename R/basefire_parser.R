.basefire_clean_input_lines <- function(lines) {
  if (!is.character(lines)) {
    .fire_stop('Base Fire input lines must be character.')
  }

  lines <- sub('^\ufeff', '', lines)
  lines <- sub('\r$', '', lines)

  # LANDIS-II parameter files commonly use lines beginning with >> as comments
  # and may contain << trailing explanatory text. We also tolerate # comments
  # for synthetic/testing files.
  lines[grepl('^\\s*>>', lines)] <- ''
  lines <- sub('<<.*$', '', lines)
  lines <- sub('#.*$', '', lines)
  lines <- trimws(lines)

  lines[nzchar(lines)]
}

.basefire_tokenize_line <- function(line) {
  scan(
    text = line,
    what = character(),
    quiet = TRUE,
    quote = c('"', "'"),
    comment.char = ''
  )
}

.basefire_numeric_token <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1L || is.na(value) || !is.finite(value)) {
    .fire_stop('Invalid numeric value for ', name, ': ', x, '.')
  }
  value
}

.basefire_integer_token <- function(x, name) {
  value <- .basefire_numeric_token(x, name)
  out <- .basefire_integer_like(value, name)
  if (length(out) != 1L) {
    .fire_stop('Invalid integer value for ', name, '.')
  }
  out
}

.basefire_percentage_token <- function(x, name) {
  if (grepl('%$', x)) {
    value <- .basefire_numeric_token(sub('%$', '', x), name) / 100
  } else {
    value <- .basefire_numeric_token(x, name)
  }
  value
}

.basefire_resolve_path <- function(path, base_dir) {
  .assert_scalar_character(path, 'path')

  is_absolute <- grepl('^[A-Za-z]:[/\\\\]', path) ||
    grepl('^[/\\\\]{2}', path) ||
    grepl('^/', path)

  if (is_absolute) path else file.path(base_dir, path)
}

.basefire_read_species_csv <- function(path) {
  if (!file.exists(path)) {
    .fire_stop('Base Fire species CSV not found: ', path)
  }

  x <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required <- c('SpeciesCode', 'FireTolerance')
  missing <- setdiff(required, names(x))

  if (length(missing)) {
    .fire_stop(
      'Base Fire species CSV is missing column(s): ',
      paste(missing, collapse = ', '),
      '.'
    )
  }

  .basefire_normalize_species(data.frame(
    species_code = as.character(x[['SpeciesCode']]),
    fire_tolerance = suppressWarnings(as.numeric(x[['FireTolerance']])),
    stringsAsFactors = FALSE
  ))
}

.basefire_read_parameters <- function(path, known_species = NULL,
                                      base_dir = dirname(path), project = NULL) {
  if (!file.exists(path)) {
    .fire_stop('Base Fire parameter file not found: ', path)
  }

  source_file <- normalizePath(path, winslash = '/', mustWork = TRUE)
  lines <- .basefire_clean_input_lines(
    readLines(path, warn = FALSE, encoding = 'UTF-8')
  )

  if (!length(lines)) {
    .fire_stop('Base Fire parameter file is empty.')
  }

  tokens <- lapply(lines, .basefire_tokenize_line)
  i <- 1L

  current_first <- function() {
    if (i > length(tokens) || !length(tokens[[i]])) return(NA_character_)
    tokens[[i]][1L]
  }

  read_named_value <- function(key) {
    if (i > length(tokens)) {
      .fire_stop('Expected ', key, ' but reached end of Base Fire input.')
    }

    tok <- tokens[[i]]
    if (!length(tok) || !identical(tok[1L], key)) {
      .fire_stop(
        'Expected Base Fire parameter ',
        key,
        ' at logical input line ',
        i,
        '.'
      )
    }

    if (length(tok) >= 2L) {
      value <- paste(tok[-1L], collapse = ' ')
      i <<- i + 1L
      return(value)
    }

    i <<- i + 1L
    if (i > length(tokens) || !length(tokens[[i]])) {
      .fire_stop('Missing value for Base Fire parameter ', key, '.')
    }

    value <- paste(tokens[[i]], collapse = ' ')
    i <<- i + 1L
    value
  }

  if (!identical(current_first(), 'LandisData')) {
    .fire_stop('Base Fire input must begin with LandisData.')
  }
  landis_data <- read_named_value('LandisData')

  timestep <- .basefire_integer_token(
    read_named_value('Timestep'),
    'Timestep'
  )

  species_csv_file <- NULL
  if (identical(current_first(), 'Species_CSV_File')) {
    species_csv_file <- read_named_value('Species_CSV_File')
    species_path <- .basefire_resolve_path(species_csv_file, base_dir)
    species <- .basefire_read_species_csv(species_path)
  } else {
    if (!inherits(project, 'landisR_project')) {
      .fire_stop(
        'Supply project = <landisR_project> to read a LANDIS-II Base Fire ',
        'file without Species_CSV_File. Fire tolerance comes from the project species table.'
      )
    }
    project_species <- project$canonical$species
    .basefire_require_columns(
      project_species, c('species', 'fire_tolerance'), 'project species table'
    )
    species <- .basefire_normalize_species(data.frame(
      species_code = project_species$species,
      fire_tolerance = project_species$fire_tolerance
    ))
  }

  if (!is.null(known_species)) {
    unknown <- setdiff(species$species_code, as.character(known_species))
    if (length(unknown)) {
      .fire_stop(
        'Base Fire species table contains species not present in the supplied ',
        'landisR species set: ', paste(unknown, collapse = ', '), '.'
      )
    }
  }

  region_rows <- list()
  while (i <= length(tokens) &&
         !identical(current_first(), 'InitialFireRegionsMap')) {
    tok <- tokens[[i]]
    if (length(tok) != 7L) {
      .fire_stop(
        'A Base Fire region row must contain 7 fields at logical input line ',
        i,
        '.'
      )
    }

    region_rows[[length(region_rows) + 1L]] <- data.frame(
      name = tok[1L],
      map_code = .basefire_integer_token(tok[2L], 'Map Code'),
      mean_size = .basefire_numeric_token(tok[3L], 'Mean Fire Size'),
      min_size = .basefire_numeric_token(tok[4L], 'Minimum Fire Size'),
      max_size = .basefire_numeric_token(tok[5L], 'Maximum Fire Size'),
      ignition_probability = .basefire_numeric_token(
        tok[6L],
        'Ignition Probability'
      ),
      fire_spread_age = .basefire_integer_token(
        tok[7L],
        'Fire Spread Age'
      ),
      stringsAsFactors = FALSE
    )

    i <- i + 1L
  }

  if (!length(region_rows)) {
    .fire_stop('No Base Fire regions were defined.')
  }

  regions <- do.call(rbind, region_rows)
  initial_fire_regions <- read_named_value('InitialFireRegionsMap')

  dynamic_regions <- .basefire_empty_dynamic_regions()

  if (i <= length(tokens) &&
      identical(current_first(), 'DynamicFireRegionTable')) {
    i <- i + 1L
    rows <- list()

    while (i <= length(tokens) &&
           !identical(current_first(), 'FuelCurveTable')) {
      tok <- tokens[[i]]
      if (length(tok) != 2L) {
        .fire_stop(
          'A DynamicFireRegionTable row must contain year and map name.'
        )
      }

      rows[[length(rows) + 1L]] <- data.frame(
        year = .basefire_integer_token(tok[1L], 'Dynamic fire-region year'),
        map_name = tok[2L],
        stringsAsFactors = FALSE
      )
      i <- i + 1L
    }

    if (length(rows)) {
      dynamic_regions <- do.call(rbind, rows)
    }
  }

  if (!identical(current_first(), 'FuelCurveTable')) {
    .fire_stop('Expected FuelCurveTable in Base Fire input.')
  }
  i <- i + 1L

  fuel_rows <- list()
  while (i <= length(tokens) &&
         !identical(current_first(), 'WindCurveTable')) {
    tok <- tokens[[i]]
    if (length(tok) != 6L) {
      .fire_stop('A FuelCurveTable row must contain region + 5 ages.')
    }

    fuel_rows[[length(fuel_rows) + 1L]] <- data.frame(
      region_name = tok[1L],
      severity1 = .basefire_integer_token(tok[2L], 'Fuel severity1'),
      severity2 = .basefire_integer_token(tok[3L], 'Fuel severity2'),
      severity3 = .basefire_integer_token(tok[4L], 'Fuel severity3'),
      severity4 = .basefire_integer_token(tok[5L], 'Fuel severity4'),
      severity5 = .basefire_integer_token(tok[6L], 'Fuel severity5'),
      stringsAsFactors = FALSE
    )
    i <- i + 1L
  }

  fuel_curves <- if (length(fuel_rows)) {
    do.call(rbind, fuel_rows)
  } else {
    .basefire_empty_curves()
  }

  if (!identical(current_first(), 'WindCurveTable')) {
    .fire_stop('Expected WindCurveTable in Base Fire input.')
  }
  i <- i + 1L

  wind_rows <- list()
  while (i <= length(tokens) &&
         !identical(current_first(), 'FireDamageTable')) {
    tok <- tokens[[i]]
    if (length(tok) != 6L) {
      .fire_stop('A WindCurveTable row must contain region + 5 ages.')
    }

    # InputParameterParser reads wind values in file order S5, S4, S3, S2, S1.
    wind_rows[[length(wind_rows) + 1L]] <- data.frame(
      region_name = tok[1L],
      severity1 = .basefire_integer_token(tok[6L], 'Wind severity1'),
      severity2 = .basefire_integer_token(tok[5L], 'Wind severity2'),
      severity3 = .basefire_integer_token(tok[4L], 'Wind severity3'),
      severity4 = .basefire_integer_token(tok[3L], 'Wind severity4'),
      severity5 = .basefire_integer_token(tok[2L], 'Wind severity5'),
      stringsAsFactors = FALSE
    )
    i <- i + 1L
  }

  wind_curves <- if (length(wind_rows)) {
    do.call(rbind, wind_rows)
  } else {
    .basefire_empty_curves()
  }

  if (!identical(current_first(), 'FireDamageTable')) {
    .fire_stop('Expected FireDamageTable in Base Fire input.')
  }
  i <- i + 1L

  damage_rows <- list()
  while (i <= length(tokens) &&
         !identical(current_first(), 'MapNames')) {
    tok <- tokens[[i]]
    if (length(tok) != 2L) {
      .fire_stop('A FireDamageTable row must contain 2 fields.')
    }

    damage_rows[[length(damage_rows) + 1L]] <- data.frame(
      max_age = .basefire_percentage_token(tok[1L], 'Max Survival Age'),
      severity_tolerance_diff = .basefire_integer_token(
        tok[2L],
        'Severity Tolerance Diff'
      ),
      stringsAsFactors = FALSE
    )
    i <- i + 1L
  }

  if (!length(damage_rows)) {
    .fire_stop('No Base Fire damage classes were defined.')
  }
  damage <- do.call(rbind, damage_rows)

  map_names_template <- read_named_value('MapNames')
  log_file <- read_named_value('LogFile')
  summary_log_file <- read_named_value('SummaryLogFile')

  if (i <= length(tokens)) {
    .fire_stop('Unexpected data after SummaryLogFile in Base Fire input.')
  }

  x <- .basefire_parameters(
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

  attr(x, 'landis_data') <- landis_data
  x
}
