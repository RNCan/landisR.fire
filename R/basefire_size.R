.basefire_default_exponential_draw <- function(lambda) {
  if (!is.numeric(lambda) ||
      length(lambda) != 1L ||
      is.na(lambda) ||
      !is.finite(lambda) ||
      lambda <= 0) {
    .fire_stop('lambda must be one finite numeric value > 0.')
  }

  stats::rexp(1L, rate = lambda)
}

.basefire_prepare_size <- function(
    landscape,
    fire_state,
    topology,
    context,
    model) {

  cache <- context$basefire_cache

  if (is.null(cache) ||
      !inherits(cache, 'landisR_fire_basefire_parameter_cache')) {
    .fire_stop(
      'Base Fire size preparation requires context$basefire_cache.'
    )
  }

  cell_area <- context$cell_area

  if (!is.numeric(cell_area) ||
      length(cell_area) != 1L ||
      is.na(cell_area) ||
      !is.finite(cell_area) ||
      cell_area <= 0) {
    .fire_stop(
      'Base Fire size preparation requires context$cell_area > 0.'
    )
  }

  x <- list(
    region_index = cache$region_index,
    basefire_index = cache$basefire_index,
    region_name = cache$region_name,
    mean_size = as.numeric(cache$mean_size),
    min_size = as.numeric(cache$min_size),
    max_size = as.numeric(cache$max_size),
    cell_area = as.numeric(cell_area)
  )

  class(x) <- 'landisR_fire_basefire_size_cache'
  .basefire_validate_size_cache(x)
  x
}

.basefire_validate_size_cache <- function(x) {
  if (!inherits(x, 'landisR_fire_basefire_size_cache')) {
    .fire_stop(
      'prepared must inherit from landisR_fire_basefire_size_cache.'
    )
  }

  n <- length(x$region_index)

  fields <- c(
    'basefire_index',
    'region_name',
    'mean_size',
    'min_size',
    'max_size'
  )

  for (nm in fields) {
    if (length(x[[nm]]) != n) {
      .fire_stop('Invalid Base Fire size cache field: ', nm, '.')
    }
  }

  if (!is.numeric(x$cell_area) ||
      length(x$cell_area) != 1L ||
      is.na(x$cell_area) ||
      !is.finite(x$cell_area) ||
      x$cell_area <= 0) {
    .fire_stop('Invalid Base Fire size cache cell_area.')
  }

  invisible(TRUE)
}

.basefire_compute_size <- function(
    region_index,
    prepared,
    exponential_draw = .basefire_default_exponential_draw,
    diagnostics = FALSE) {

  .basefire_validate_size_cache(prepared)
  .assert_function(exponential_draw, 'exponential_draw')

  region_index <- as.integer(region_index)

  if (length(region_index) != 1L ||
      is.na(region_index) ||
      region_index < 1L ||
      region_index > length(prepared$region_index)) {
    .fire_stop('region_index must identify one prepared Base Fire region.')
  }

  mean_size <- prepared$mean_size[region_index]
  min_size <- prepared$min_size[region_index]
  max_size <- prepared$max_size[region_index]
  cell_area <- prepared$cell_area

  result <- list(
    region_index = region_index,
    basefire_index = prepared$basefire_index[region_index],
    region_name = prepared$region_name[region_index],
    mean_size = mean_size,
    min_size = min_size,
    max_size = max_size,
    cell_area = cell_area,
    lambda = NA_real_,
    size_generated = NA_real_,
    bounded_size = NA_real_,
    size_sites = 0L,
    draws = 0L,
    bound = 'mean_nonpositive'
  )

  # FireEvent.ComputeSize returns immediately and consumes no exponential
  # random draw when MeanSize <= 0.
  if (mean_size <= 0) {
    class(result) <- 'landisR_fire_basefire_size_result'
    return(if (isTRUE(diagnostics)) result else 0L)
  }

  lambda <- 1.0 / mean_size
  size_generated <- exponential_draw(lambda)

  if (!is.numeric(size_generated) ||
      length(size_generated) != 1L ||
      is.na(size_generated) ||
      !is.finite(size_generated) ||
      size_generated < 0) {
    .fire_stop(
      'exponential_draw(lambda) must return one finite value >= 0.'
    )
  }

  # Exact C# branch structure:
  #   if generated < MinSize -> use MinSize
  #   else if generated > MaxSize -> use MaxSize
  #   else -> use generated
  #
  # The reference parser does not require MinSize <= MeanSize <= MaxSize,
  # so we intentionally do not normalize or reorder these values here.
  if (size_generated < min_size) {
    bounded_size <- min_size
    bound <- 'min'
  } else if (size_generated > max_size) {
    bounded_size <- max_size
    bound <- 'max'
  } else {
    bounded_size <- size_generated
    bound <- 'none'
  }

  # C# casts positive doubles to int by truncating toward zero.
  # All valid Base Fire sizes are non-negative, so floor() is equivalent.
  size_sites_double <- floor(bounded_size / cell_area)

  if (!is.finite(size_sites_double) ||
      size_sites_double > .Machine$integer.max) {
    .fire_stop(
      'Computed Base Fire size exceeds R integer capacity.'
    )
  }

  size_sites <- as.integer(size_sites_double)

  result$lambda <- lambda
  result$size_generated <- as.numeric(size_generated)
  result$bounded_size <- as.numeric(bounded_size)
  result$size_sites <- size_sites
  result$draws <- 1L
  result$bound <- bound

  class(result) <- 'landisR_fire_basefire_size_result'

  if (isTRUE(diagnostics)) result else size_sites
}

.basefire_size <- .wrap_component(
  fun = function(
      region_index,
      prepared,
      exponential_draw = .basefire_default_exponential_draw,
      diagnostics = FALSE,
      ...) {

    .basefire_compute_size(
      region_index = region_index,
      prepared = prepared,
      exponential_draw = exponential_draw,
      diagnostics = diagnostics
    )
  },
  class_name = 'size',
  prepare = .basefire_prepare_size,
  metadata = list(
    profile = 'basefire',
    reference = 'FireEvent.ComputeSize'
  )
)
