#' Build a compact regular-grid topology for fire spread
#'
#' The topology stores internal integer indices rather than repeatedly joining
#' or searching by cell identifier inside fire-spread loops. Missing/inactive
#' neighbors are represented by 0L. The first and second directional rings are
#' precomputed for eight directions.
#'
#' @param landscape data.frame containing cell identifier, row and column.
#' @param cell_id Name of the cell identifier column.
#' @param row Name of the row column.
#' @param col Name of the column column.
#' @param active Optional active-mask column. If absent, every row is active.
#' @export
fire_grid_topology <- function(landscape, cell_id = 'cell_id', row = 'row',
                               col = 'col', active = 'active') {
  if (!is.data.frame(landscape)) {
    .fire_stop('landscape must be a data.frame.')
  }

  required <- c(cell_id, row, col)
  missing <- setdiff(required, names(landscape))
  if (length(missing)) {
    .fire_stop(
      'landscape is missing required columns: ',
      paste(missing, collapse = ', '),
      '.'
    )
  }

  keep <- rep(TRUE, nrow(landscape))

  if (active %in% names(landscape)) {
    keep <- as.logical(landscape[[active]])
    if (anyNA(keep)) {
      .fire_stop('active mask cannot contain NA values.')
    }
  }

  x <- landscape[keep, required, drop = FALSE]

  ids <- as.integer(x[[cell_id]])
  rows <- as.integer(x[[row]])
  cols <- as.integer(x[[col]])

  if (!length(ids) || anyNA(ids) || anyDuplicated(ids)) {
    .fire_stop('active cell identifiers must be unique non-missing integers.')
  }
  if (anyNA(rows) || anyNA(cols)) {
    .fire_stop('row and col must be finite integer-like values.')
  }
  if (anyDuplicated(data.frame(row = rows, col = cols))) {
    .fire_stop('active row/col coordinates must be unique.')
  }

  row_min <- min(rows)
  col_min <- min(cols)
  row_internal <- rows - row_min + 1L
  col_internal <- cols - col_min + 1L

  n_row <- max(row_internal)
  n_col <- max(col_internal)

  lookup <- matrix(
    0L,
    nrow = n_row,
    ncol = n_col
  )

  lookup[cbind(row_internal, col_internal)] <- seq_along(ids)

  directions <- data.frame(
    direction = c('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'),
    drow = c(-1L, -1L, 0L, 1L, 1L, 1L, 0L, -1L),
    dcol = c(0L, 1L, 1L, 1L, 0L, -1L, -1L, -1L),
    stringsAsFactors = FALSE
  )

  build_ring <- function(multiplier) {
    ans <- matrix(
      0L,
      nrow = length(ids),
      ncol = nrow(directions),
      dimnames = list(NULL, directions$direction)
    )

    for (j in seq_len(nrow(directions))) {
      rr <- row_internal + multiplier * directions$drow[j]
      cc <- col_internal + multiplier * directions$dcol[j]

      valid <- rr >= 1L & rr <= n_row & cc >= 1L & cc <= n_col

      if (any(valid)) {
        ans[valid, j] <- lookup[cbind(rr[valid], cc[valid])]
      }
    }

    ans
  }

  topology <- list(
    cell_id = ids,
    row = rows,
    col = cols,
    internal_index = seq_along(ids),
    row_min = row_min,
    col_min = col_min,
    n_row = n_row,
    n_col = n_col,
    lookup = lookup,
    directions = directions,
    neighbors = build_ring(1L),
    second_neighbors = build_ring(2L)
  )

  class(topology) <- 'landisR_fire_topology'
  validate_fire_topology(topology)
  topology
}

#' Validate a fire-grid topology
#' @param x Object to validate.
#' @export
validate_fire_topology <- function(x) {
  if (!inherits(x, 'landisR_fire_topology')) {
    .fire_stop('x must inherit from landisR_fire_topology.')
  }

  n <- length(x$cell_id)

  if (!n || anyDuplicated(x$cell_id)) {
    .fire_stop('topology cell identifiers must be unique.')
  }
  if (!is.matrix(x$neighbors) || nrow(x$neighbors) != n ||
      ncol(x$neighbors) != 8L) {
    .fire_stop('topology neighbors must be an n-cell by 8 integer matrix.')
  }
  if (!is.matrix(x$second_neighbors) || nrow(x$second_neighbors) != n ||
      ncol(x$second_neighbors) != 8L) {
    .fire_stop('topology second_neighbors must be an n-cell by 8 matrix.')
  }
  if (any(x$neighbors < 0L) || any(x$neighbors > n) ||
      any(x$second_neighbors < 0L) || any(x$second_neighbors > n)) {
    .fire_stop('topology neighbor indices must be 0 or valid internal indices.')
  }

  invisible(TRUE)
}

#' Prepare reusable fire-model caches
#'
#' Component preparation is deliberately separated from event execution. A
#' component may attach a `prepare` function through one of the custom wrapper
#' constructors. fire_prepare() calls each such function once and stores its
#' result in `prepared[[component_name]]`.
#'
#' A preparation function must accept the formal arguments
#' `landscape`, `fire_state`, `topology`, `context`, and `model` (or `...`).
#' Runtime component functions can later receive the corresponding prepared
#' cache directly, avoiding repeated joins or extraction inside cell loops.
#'
#' @param model A landisR_fire_model.
#' @param landscape Fire-model landscape view or other model-specific object.
#' @param fire_state Optional external fire state.
#' @param context Optional named/list context passed to component preparers.
#' @param topology Optional prebuilt landisR_fire_topology. If NULL and
#'   landscape contains cell_id, row and col, topology is built automatically.
#' @export
fire_prepare <- function(model, landscape, fire_state = NULL,
                         context = list(), topology = NULL) {
  validate_fire_model(model)
  .assert_named_list(context, 'context')

  if (is.null(topology) && is.data.frame(landscape) &&
      all(c('cell_id', 'row', 'col') %in% names(landscape))) {
    topology <- fire_grid_topology(landscape)
  }

  if (!is.null(topology)) {
    validate_fire_topology(topology)
  }

  prepared <- setNames(
    vector('list', length(.fire_component_names)),
    .fire_component_names
  )

  for (nm in .fire_component_names) {
    prepare_fun <- attr(
      model[[nm]],
      'fire_prepare',
      exact = TRUE
    )

    if (!is.null(prepare_fun)) {
      prepared[[nm]] <- prepare_fun(
        landscape = landscape,
        fire_state = fire_state,
        topology = topology,
        context = context,
        model = model
      )
    }
  }

  x <- model
  x$prepared <- prepared
  x$topology <- topology
  x$prepare_context <- context
  x$prepared_n_cells <- if (is.null(topology)) {
    NA_integer_
  } else {
    length(topology$cell_id)
  }

  class(x) <- c(
    'landisR_fire_prepared_model',
    'landisR_fire_model'
  )

  validate_fire_prepared_model(x)
  x
}

#' Validate a prepared fire model
#' @param x Object to validate.
#' @export
validate_fire_prepared_model <- function(x) {
  if (!inherits(x, 'landisR_fire_prepared_model')) {
    .fire_stop('x must inherit from landisR_fire_prepared_model.')
  }

  validate_fire_model(x)

  if (!is.list(x$prepared) ||
      !identical(names(x$prepared), .fire_component_names)) {
    .fire_stop('prepared component cache is invalid.')
  }

  if (!is.null(x$topology)) {
    validate_fire_topology(x$topology)
  }

  invisible(TRUE)
}

#' @export
print.landisR_fire_topology <- function(x, ...) {
  cat('<landisR_fire_topology>\n')
  cat('  active cells : ', length(x$cell_id), '\n', sep = '')
  cat('  grid         : ', x$n_row, ' x ', x$n_col, '\n', sep = '')
  cat('  directions   : 8 first-ring + 8 second-ring lookups\n')
  invisible(x)
}

#' @export
print.landisR_fire_prepared_model <- function(x, ...) {
  cat('<landisR_fire_prepared_model>\n')
  cat('  name       : ', x$name, '\n', sep = '')
  cat('  cells      : ', x$prepared_n_cells, '\n', sep = '')

  cached <- vapply(
    x$prepared,
    function(z) !is.null(z),
    logical(1)
  )

  cat(
    '  prepared   : ',
    if (any(cached)) paste(names(cached)[cached], collapse = ', ') else 'none',
    '\n',
    sep = ''
  )
  invisible(x)
}
