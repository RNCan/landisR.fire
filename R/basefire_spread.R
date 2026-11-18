.basefire_wind_probabilities <- function(wind_speed) {
  if (!is.numeric(wind_speed) ||
      length(wind_speed) != 1L ||
      is.na(wind_speed) ||
      !is.finite(wind_speed) ||
      wind_speed < 0 ||
      wind_speed > 1) {
    .fire_stop('wind_speed must be one finite value between 0 and 1.')
  }

  base <- (4.0 - wind_speed) / 8.0

  c(
    base * (1.0 + wind_speed),
    base * (1.0 + wind_speed),
    base,
    base * (1.0 - wind_speed),
    base * (1.0 - wind_speed),
    base * (1.0 - wind_speed),
    base,
    base * (1.0 + wind_speed)
  )
}

.basefire_site_directions <- function() {
  data.frame(
    direction = c('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'),
    drow = c(-1L, -1L, 0L, 1L, 1L, 1L, 0L, -1L),
    dcol = c(0L, 1L, 1L, 1L, 0L, -1L, -1L, -1L),
    stringsAsFactors = FALSE
  )
}

.basefire_shared_geometry <- function(landscape) {
  if (anyNA(landscape$active))
    .fire_stop("The Base Fire active mask cannot contain NA.")
  if (anyNA(landscape$cell_id) || anyDuplicated(landscape$cell_id))
    .fire_stop("Base Fire landscape cell identifiers must be unique and non-missing.")

  active_to_full <- which(landscape$active)
  if (!length(active_to_full))
    .fire_stop("Base Fire site topology cannot be empty.")
  active <- structure(list(
    cell_id = landscape$cell_id[active_to_full],
    row = landscape$row[active_to_full],
    col = landscape$col[active_to_full]
  ), class = "landisR_basefire_active_cells")

  full_to_active <- integer(nrow(landscape))
  full_to_active[active_to_full] <- seq_along(active_to_full)
  directions <- .basefire_site_directions()
  rings <- basefire_site_neighbors(
    landscape$row, landscape$col, active_to_full,
    directions$drow, directions$dcol, directions$direction
  )
  sites <- structure(list(
    active_cell_id = active$cell_id,
    full_cell_id = landscape$cell_id,
    full_active = landscape$active,
    full_to_active = full_to_active,
    active_to_full = active_to_full,
    directions = directions,
    scan_order_v7 = order(active$row, active$col, method = "radix"),
    neighbors_full = rings[[1L]],
    second_neighbors_full = rings[[2L]]
  ), class = "landisR_fire_basefire_site_topology")
  list(active = active, sites = sites)
}

.basefire_build_site_topology <- function(
    landscape,
    active_topology,
    full_landscape = NULL) {

  if (is.null(full_landscape)) {
    full_landscape <- landscape
  }

  if (!is.data.frame(full_landscape)) {
    .fire_stop('full Base Fire landscape must be a data.frame.')
  }

  required <- c('cell_id', 'row', 'col')
  missing <- setdiff(required, names(full_landscape))

  if (length(missing)) {
    .fire_stop(
      'full Base Fire landscape is missing columns: ',
      paste(missing, collapse = ', '),
      '.'
    )
  }

  validate_fire_topology(active_topology)

  full_cell_id <- .basefire_integer_like(
    full_landscape$cell_id,
    'full_landscape$cell_id'
  )
  full_row <- .basefire_integer_like(
    full_landscape$row,
    'full_landscape$row'
  )
  full_col <- .basefire_integer_like(
    full_landscape$col,
    'full_landscape$col'
  )

  if (anyDuplicated(full_cell_id)) {
    .fire_stop('full Base Fire landscape cell_id values must be unique.')
  }

  if (anyDuplicated(data.frame(row = full_row, col = full_col))) {
    .fire_stop('full Base Fire landscape row/col coordinates must be unique.')
  }

  if ('active' %in% names(full_landscape)) {
    full_active <- as.logical(full_landscape$active)

    if (anyNA(full_active)) {
      .fire_stop('full Base Fire landscape active mask cannot contain NA.')
    }
  } else {
    full_active <- rep(TRUE, length(full_cell_id))
  }

  active_to_full <- match(active_topology$cell_id, full_cell_id)

  if (anyNA(active_to_full)) {
    .fire_stop(
      'Every active topology cell must occur in the full Base Fire landscape.'
    )
  }

  if (any(!full_active[active_to_full])) {
    .fire_stop(
      'An active topology cell is marked inactive in the full landscape.'
    )
  }

  full_to_active <- match(full_cell_id, active_topology$cell_id)
  full_to_active[is.na(full_to_active)] <- 0L
  full_to_active <- as.integer(full_to_active)

  inconsistent <- full_active != (full_to_active > 0L)

  if (any(inconsistent)) {
    .fire_stop(
      'The active mask of the full landscape does not match the active ',
      'fire topology.'
    )
  }

  row_min <- min(full_row)
  col_min <- min(full_col)

  row_internal <- full_row - row_min + 1L
  col_internal <- full_col - col_min + 1L

  n_row <- max(row_internal)
  n_col <- max(col_internal)

  lookup <- matrix(
    0L,
    nrow = n_row,
    ncol = n_col
  )

  lookup[cbind(row_internal, col_internal)] <- seq_along(full_cell_id)

  source_row <- row_internal[active_to_full]
  source_col <- col_internal[active_to_full]

  directions <- .basefire_site_directions()

  build_ring <- function(multiplier) {
    out <- matrix(
      0L,
      nrow = length(active_topology$cell_id),
      ncol = 8L,
      dimnames = list(NULL, directions$direction)
    )

    for (j in seq_len(8L)) {
      rr <- source_row + multiplier * directions$drow[j]
      cc <- source_col + multiplier * directions$dcol[j]

      valid <- rr >= 1L & rr <= n_row & cc >= 1L & cc <= n_col

      if (any(valid)) {
        out[valid, j] <- lookup[cbind(rr[valid], cc[valid])]
      }
    }

    out
  }

  x <- list(
    active_cell_id = active_topology$cell_id,
    full_cell_id = full_cell_id,
    full_active = full_active,
    full_to_active = full_to_active,
    active_to_full = as.integer(active_to_full),
    directions = directions,
    scan_order_v7 = order(source_row, source_col, method = "radix"),
    neighbors_full = build_ring(1L),
    second_neighbors_full = build_ring(2L)
  )

  class(x) <- 'landisR_fire_basefire_site_topology'
  .basefire_validate_site_topology(x)
  x
}

.basefire_validate_site_topology <- function(x) {
  if (!inherits(x, 'landisR_fire_basefire_site_topology')) {
    .fire_stop(
      'x must inherit from landisR_fire_basefire_site_topology.'
    )
  }

  n_active <- length(x$active_cell_id)
  n_full <- length(x$full_cell_id)

  if (!n_active || !n_full) {
    .fire_stop('Base Fire site topology cannot be empty.')
  }

  if (!is.matrix(x$neighbors_full) ||
      nrow(x$neighbors_full) != n_active ||
      ncol(x$neighbors_full) != 8L) {
    .fire_stop('neighbors_full must be n-active by 8.')
  }

  if (!is.matrix(x$second_neighbors_full) ||
      nrow(x$second_neighbors_full) != n_active ||
      ncol(x$second_neighbors_full) != 8L) {
    .fire_stop('second_neighbors_full must be n-active by 8.')
  }

  if (any(x$neighbors_full < 0L | x$neighbors_full > n_full) ||
      any(x$second_neighbors_full < 0L |
          x$second_neighbors_full > n_full)) {
    .fire_stop('Base Fire full-neighbor indices are invalid.')
  }

  if (length(x$full_active) != n_full ||
      length(x$full_to_active) != n_full ||
      length(x$active_to_full) != n_active) {
    .fire_stop('Base Fire site topology mapping lengths are invalid.')
  }

  invisible(TRUE)
}

.basefire_prepare_spread <- function(
    landscape,
    fire_state,
    topology,
    context,
    model) {

  if (is.null(topology)) {
    topology <- fire_grid_topology(landscape)
  }

  .validate_fire_preparation_topology(topology)

  cache <- context$basefire_cache

  if (is.null(cache) ||
      !inherits(cache, 'landisR_fire_basefire_parameter_cache')) {
    .fire_stop(
      'Base Fire spread preparation requires context$basefire_cache.'
    )
  }

  site_topology <- context$basefire_site_topology
  if (is.null(site_topology)) {
    site_topology <- .basefire_build_site_topology(
      landscape = landscape,
      active_topology = topology,
      full_landscape = context$basefire_full_landscape
    )
  }

  region_index <- .basefire_resolve_region_index(
    landscape = landscape,
    topology = topology,
    cache = cache
  )

  x <- list(
    cell_id = topology$cell_id,
    region_index = region_index,
    n_regions = cache$n_regions,
    site_topology = site_topology
  )

  class(x) <- 'landisR_fire_basefire_spread_cache'
  .basefire_validate_spread_cache(x)
  x
}

.basefire_validate_spread_cache <- function(x) {
  if (!inherits(x, 'landisR_fire_basefire_spread_cache')) {
    .fire_stop(
      'prepared must inherit from landisR_fire_basefire_spread_cache.'
    )
  }

  .basefire_validate_site_topology(x$site_topology)

  n <- length(x$cell_id)

  if (length(x$region_index) != n ||
      any(x$region_index < 1L | x$region_index > x$n_regions)) {
    .fire_stop('Invalid Base Fire spread region index.')
  }

  invisible(TRUE)
}

.basefire_get_neighbors <- function(
    cell_index,
    wind_direction,
    wind_speed,
    prepared,
    next_double = stats::runif,
    generate_uniform = stats::runif) {

  .basefire_validate_spread_cache(prepared)
  .assert_function(next_double, 'next_double')
  .assert_function(generate_uniform, 'generate_uniform')

  cell_index <- as.integer(cell_index)
  wind_direction <- as.integer(wind_direction)

  if (length(cell_index) != 1L ||
      is.na(cell_index) ||
      cell_index < 1L ||
      cell_index > length(prepared$cell_id)) {
    .fire_stop('cell_index must identify one active Base Fire cell.')
  }

  if (length(wind_direction) != 1L ||
      is.na(wind_direction) ||
      wind_direction < 0L) {
    .fire_stop('wind_direction must be one integer >= 0.')
  }

  # Literal GetNeighbors() behavior.
  if (wind_direction > 7L) {
    wind_direction <- 7L
  }

  wind_probs <- .basefire_wind_probabilities(wind_speed)
  topo <- prepared$site_topology

  result <- vector('list', 9L)
  n_result <- 0L
  next_double_draws <- 0L
  uniform_draws <- 0L

  for (j in seq_len(8L)) {
    full_index <- topo$neighbors_full[cell_index, j]

    # C# short-circuit:
    # if neighbor == null, NextDouble() is NOT consumed.
    if (full_index == 0L) {
      next
    }

    prob_index <- ((j - 1L + wind_direction) %% 8L) + 1L
    probability <- wind_probs[prob_index]

    draw <- .basefire_uniform_one(next_double)
    next_double_draws <- next_double_draws + 1L

    if (draw < probability) {
      n_result <- n_result + 1L

      result[[n_result]] <- data.frame(
        full_index = as.integer(full_index),
        active_index = topo$full_to_active[full_index],
        cell_id = topo$full_cell_id[full_index],
        is_active = topo$full_active[full_index],
        ring = 1L,
        direction_index = j - 1L,
        direction = topo$directions$direction[j],
        selection_probability = probability,
        selection_draw = draw,
        stringsAsFactors = FALSE
      )
    }
  }

  # Ninth neighbor is two cells away in the selected wind direction.
  full_index9 <- topo$second_neighbors_full[
    cell_index,
    wind_direction + 1L
  ]

  if (full_index9 != 0L) {
    draw9 <- .basefire_uniform_one(generate_uniform)
    uniform_draws <- uniform_draws + 1L

    if (draw9 < wind_speed) {
      n_result <- n_result + 1L

      result[[n_result]] <- data.frame(
        full_index = as.integer(full_index9),
        active_index = topo$full_to_active[full_index9],
        cell_id = topo$full_cell_id[full_index9],
        is_active = topo$full_active[full_index9],
        ring = 2L,
        direction_index = wind_direction,
        direction = topo$directions$direction[wind_direction + 1L],
        selection_probability = wind_speed,
        selection_draw = draw9,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!n_result) {
    out <- data.frame(
      full_index = integer(),
      active_index = integer(),
      cell_id = integer(),
      is_active = logical(),
      ring = integer(),
      direction_index = integer(),
      direction = character(),
      selection_probability = numeric(),
      selection_draw = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, result[seq_len(n_result)])
    rownames(out) <- NULL
  }

  attr(out, 'next_double_draws') <- next_double_draws
  attr(out, 'uniform_draws') <- uniform_draws
  attr(out, 'wind_direction') <- wind_direction
  attr(out, 'wind_speed') <- wind_speed
  out
}

.basefire_prepare_burn_probability <- function(
    landscape,
    fire_state,
    topology,
    context,
    model) {

  ignition_cache <- .basefire_prepare_ignition(
    landscape = landscape,
    fire_state = fire_state,
    topology = topology,
    context = context,
    model = model
  )

  x <- list(
    cell_id = ignition_cache$cell_id,
    region_index = ignition_cache$region_index,
    fire_spread_age = ignition_cache$fire_spread_age,
    last_fire_year = ignition_cache$last_fire_year
  )

  class(x) <- 'landisR_fire_basefire_burn_probability_cache'
  .basefire_validate_burn_probability_cache(x)
  x
}

.basefire_validate_burn_probability_cache <- function(x) {
  if (!inherits(x, 'landisR_fire_basefire_burn_probability_cache')) {
    .fire_stop(
      'prepared must inherit from ',
      'landisR_fire_basefire_burn_probability_cache.'
    )
  }

  n <- length(x$cell_id)

  if (length(x$region_index) != n ||
      length(x$fire_spread_age) != n ||
      length(x$last_fire_year) != n) {
    .fire_stop('Invalid Base Fire burn-probability cache.')
  }

  invisible(TRUE)
}

.basefire_compute_burn_probability <- function(
    cell_index,
    current_time,
    prepared,
    last_fire_year = NULL) {

  .basefire_validate_burn_probability_cache(prepared)

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

  tsf <- current_time - last_fire_year[cell_index]

  .basefire_compute_init_spread_prob(
    time_since_last_fire = tsf,
    fire_spread_age = prepared$fire_spread_age[cell_index]
  )
}

.basefire_default_shuffle <- function(x) {
  if (!is.data.frame(x)) {
    .fire_stop('Base Fire shuffle input must be a data.frame.')
  }

  if (nrow(x) <= 1L) {
    return(x)
  }

  x[sample.int(nrow(x)), , drop = FALSE]
}

.basefire_call_size <- function(
    region_index,
    size_fun,
    size_prepared,
    exponential_draw) {

  .assert_function(size_fun, 'size_fun')

  value <- size_fun(
    region_index = region_index,
    prepared = size_prepared,
    exponential_draw = exponential_draw,
    diagnostics = FALSE
  )

  if (!is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value < 0 ||
      abs(value - round(value)) > sqrt(.Machine$double.eps)) {
    .fire_stop('size_fun must return one non-negative integer-like value.')
  }

  as.integer(round(value))
}

.basefire_run_spread_event <- function(
    initiation_cell_index,
    current_time,
    spread_prepared,
    size_prepared,
    burn_prepared,
    severity_fun,
    damage_fun = function(cell_index, severity, current_time) 0L,
    size_fun = .basefire_size,
    spread_fun = .basefire_spread,
    burn_probability_fun = .basefire_burn_probability,
    exponential_draw = .basefire_default_exponential_draw,
    generate_uniform = stats::runif,
    next_double = stats::runif,
    shuffle_fun = .basefire_default_shuffle,
    event_occupied = NULL,
    last_fire_year = NULL,
    wind_direction = NULL,
    wind_speed = NULL) {

  .basefire_validate_spread_cache(spread_prepared)
  .basefire_validate_size_cache(size_prepared)
  .basefire_validate_burn_probability_cache(burn_prepared)

  .assert_function(severity_fun, 'severity_fun')
  .assert_function(damage_fun, 'damage_fun')
  .assert_function(size_fun, 'size_fun')
  .assert_function(spread_fun, 'spread_fun')
  .assert_function(burn_probability_fun, 'burn_probability_fun')
  .assert_function(exponential_draw, 'exponential_draw')
  .assert_function(generate_uniform, 'generate_uniform')
  .assert_function(next_double, 'next_double')
  .assert_function(shuffle_fun, 'shuffle_fun')

  n <- length(spread_prepared$cell_id)
  n_regions <- spread_prepared$n_regions

  initiation_cell_index <- as.integer(initiation_cell_index)

  if (length(initiation_cell_index) != 1L ||
      is.na(initiation_cell_index) ||
      initiation_cell_index < 1L ||
      initiation_cell_index > n) {
    .fire_stop('initiation_cell_index is invalid.')
  }

  if (is.null(event_occupied)) {
    event_occupied <- rep(FALSE, n)
  }

  if (!is.logical(event_occupied) ||
      length(event_occupied) != n ||
      anyNA(event_occupied)) {
    .fire_stop('event_occupied must be one logical value per active cell.')
  }

  # Strict Base Fire compatibility: PlugIn.Run() calls Initiate() on every
  # active site without checking SiteVars.Event first. Spread() then assigns
  # SiteVars.Event[site] = this. Therefore an already occupied initiation
  # site is not rejected here; its event reference is effectively overwritten.
  # Neighbor candidates remain blocked by event_occupied as before.

  if (is.null(last_fire_year)) {
    last_fire_year <- burn_prepared$last_fire_year
  }

  if (!is.numeric(last_fire_year) ||
      length(last_fire_year) != n) {
    .fire_stop('last_fire_year must match active-cell count.')
  }

  if (is.null(wind_direction)) {
    wind_direction_raw <- as.integer(
      .basefire_uniform_one(generate_uniform) * 8.0
    )
  } else {
    wind_direction_raw <- as.integer(wind_direction)
  }

  if (length(wind_direction_raw) != 1L ||
      is.na(wind_direction_raw) ||
      wind_direction_raw < 0L) {
    .fire_stop('wind_direction must be one integer >= 0.')
  }

  wind_direction_used <- min(wind_direction_raw, 7L)

  if (is.null(wind_speed)) {
    wind_speed <- .basefire_uniform_one(generate_uniform)
  }

  .basefire_wind_probabilities(wind_speed)

  size_targets <- integer(n_regions)
  sites_in_event <- integer(n_regions)

  initial_region <- spread_prepared$region_index[initiation_cell_index]
  eco_index <- initial_region

  size_targets[eco_index] <- .basefire_call_size(
    region_index = eco_index,
    size_fun = size_fun,
    size_prepared = size_prepared,
    exponential_draw = exponential_draw
  )

  max_fire_region_size <- size_targets[eco_index]

  queue <- integer(n)
  queue_head <- 1L
  queue_tail <- 1L
  queue[1L] <- initiation_cell_index

  in_queue <- rep(FALSE, n)
  in_queue[initiation_cell_index] <- TRUE

  event_member <- rep(FALSE, n)
  burned <- rep(FALSE, n)
  disturbed <- rep(FALSE, n)
  cell_severity <- rep(NA_real_, n)

  checked_order <- integer(n)
  num_site_checked <- 0L

  total_sites_in_event <- 0L
  total_site_severities <- 0
  total_sites_damaged <- 0L
  cohorts_killed <- 0L

  while (
    queue_head <= queue_tail &&
    sites_in_event[eco_index] < size_targets[eco_index] &&
    total_sites_in_event < max_fire_region_size
  ) {
    num_site_checked <- num_site_checked + 1L

    site <- queue[queue_head]
    queue_head <- queue_head + 1L
    in_queue[site] <- FALSE

    checked_order[num_site_checked] <- site

    site_region <- spread_prepared$region_index[site]

    if (site_region != eco_index) {
      eco_index <- site_region

      # Literal Base Fire behavior: target sizes below 1 are treated as
      # uninitialized and ComputeSize() is called again on entry.
      if (size_targets[eco_index] < 1L) {
        size_targets[eco_index] <- .basefire_call_size(
          region_index = eco_index,
          size_fun = size_fun,
          size_prepared = size_prepared,
          exponential_draw = exponential_draw
        )

        if (size_targets[eco_index] > max_fire_region_size) {
          max_fire_region_size <- size_targets[eco_index]
        }
      }
    }

    # SiteVars.Event[site] = this occurs before severity is calculated.
    event_occupied[site] <- TRUE
    event_member[site] <- TRUE

    severity <- severity_fun(
      cell_index = site,
      current_time = current_time
    )

    if (!is.numeric(severity) ||
        length(severity) != 1L ||
        is.na(severity) ||
        !is.finite(severity) ||
        severity < 0) {
      .fire_stop(
        'severity_fun must return one finite non-negative value.'
      )
    }

    cell_severity[site] <- severity

    if (severity > 0) {
      sites_in_event[eco_index] <- sites_in_event[eco_index] + 1L
      total_sites_in_event <- total_sites_in_event + 1L
      total_site_severities <- total_site_severities + severity

      burned[site] <- TRUE
      disturbed[site] <- TRUE
      last_fire_year[site] <- current_time

      killed <- damage_fun(
        cell_index = site,
        severity = severity,
        current_time = current_time
      )

      if (!is.numeric(killed) ||
          length(killed) != 1L ||
          is.na(killed) ||
          !is.finite(killed) ||
          killed < 0 ||
          abs(killed - round(killed)) > sqrt(.Machine$double.eps)) {
        .fire_stop(
          'damage_fun must return one non-negative integer-like value.'
        )
      }

      killed <- as.integer(round(killed))
      cohorts_killed <- cohorts_killed + killed

      if (killed > 0L) {
        total_sites_damaged <- total_sites_damaged + 1L
      }

      neighbors <- spread_fun(
        cell_index = site,
        wind_direction = wind_direction_raw,
        wind_speed = wind_speed,
        prepared = spread_prepared,
        next_double = next_double,
        generate_uniform = generate_uniform
      )

      if (!is.data.frame(neighbors)) {
        .fire_stop('spread_fun must return a candidate data.frame.')
      }

      if (nrow(neighbors) > 0L) {
        neighbors <- shuffle_fun(neighbors)

        if (!is.data.frame(neighbors) ||
            !all(
              c('active_index', 'is_active') %in% names(neighbors)
            )) {
          .fire_stop(
            'shuffle_fun must return the candidate data.frame.'
          )
        }

        for (k in seq_len(nrow(neighbors))) {
          if (!isTRUE(neighbors$is_active[k])) {
            next
          }

          neighbor <- as.integer(neighbors$active_index[k])

          if (neighbor < 1L || neighbor > n) {
            .fire_stop('Invalid active neighbor index after spread.')
          }

          if (event_occupied[neighbor]) {
            next
          }

          if (in_queue[neighbor]) {
            next
          }

          burn_probability <- burn_probability_fun(
            cell_index = neighbor,
            current_time = current_time,
            prepared = burn_prepared,
            last_fire_year = last_fire_year
          )

          if (!is.numeric(burn_probability) ||
              length(burn_probability) != 1L) {
            .fire_stop(
              'burn_probability_fun must return one numeric value.'
            )
          }

          draw <- .basefire_uniform_one(generate_uniform)

          if (isTRUE(draw <= burn_probability)) {
            queue_tail <- queue_tail + 1L
            queue[queue_tail] <- neighbor
            in_queue[neighbor] <- TRUE
          }
        }
      }
    }
  }

  event_severity <- if (total_sites_damaged == 0L) {
    0
  } else {
    total_site_severities / total_sites_in_event
  }

  member_index <- which(event_member)

  cells <- data.frame(
    cell_index = member_index,
    cell_id = spread_prepared$cell_id[member_index],
    region_index = spread_prepared$region_index[member_index],
    severity = cell_severity[member_index],
    burned = burned[member_index],
    disturbed = disturbed[member_index],
    stringsAsFactors = FALSE
  )

  out <- list(
    initiation_cell_index = initiation_cell_index,
    initiation_cell_id = spread_prepared$cell_id[initiation_cell_index],
    wind_direction_raw = wind_direction_raw,
    wind_direction = wind_direction_used,
    wind_speed = wind_speed,
    size_targets = size_targets,
    sites_in_event = sites_in_event,
    max_fire_region_size = max_fire_region_size,
    num_site_checked = num_site_checked,
    checked_order = checked_order[seq_len(num_site_checked)],
    total_sites_in_event = total_sites_in_event,
    total_sites_damaged = total_sites_damaged,
    cohorts_killed = cohorts_killed,
    severity = event_severity,
    cells = cells,
    event_member = event_member,
    event_occupied = event_occupied,
    last_fire_year = last_fire_year
  )

  class(out) <- 'landisR_fire_basefire_spread_event'
  out
}

.basefire_spread <- .wrap_component(
  fun = function(
      cell_index,
      wind_direction,
      wind_speed,
      prepared,
      next_double = stats::runif,
      generate_uniform = stats::runif,
      ...) {

    .basefire_get_neighbors(
      cell_index = cell_index,
      wind_direction = wind_direction,
      wind_speed = wind_speed,
      prepared = prepared,
      next_double = next_double,
      generate_uniform = generate_uniform
    )
  },
  class_name = 'spread',
  prepare = .basefire_prepare_spread,
  metadata = list(
    profile = 'basefire',
    reference = 'FireEvent.GetNeighbors'
  )
)

.basefire_burn_probability <- .wrap_component(
  fun = function(
      cell_index,
      current_time,
      prepared,
      last_fire_year = NULL,
      ...) {

    .basefire_compute_burn_probability(
      cell_index = cell_index,
      current_time = current_time,
      prepared = prepared,
      last_fire_year = last_fire_year
    )
  },
  class_name = 'burn_probability',
  prepare = .basefire_prepare_burn_probability,
  metadata = list(
    profile = 'basefire',
    reference = 'FireEvent.ComputeFireInitSpreadProb'
  )
)
