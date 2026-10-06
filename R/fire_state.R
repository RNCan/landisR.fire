#' Construct external fire-model state
#'
#' @param cell_id Active cell identifiers.
#' @param last_fire_year Last fire year per cell, default NA.
#' @param last_severity Last fire severity per cell, default NA.
#' @param extra Additional model-specific state.
#' @export
fire_state <- function(cell_id, last_fire_year = NA_integer_,
                       last_severity = NA_real_, extra = list()) {
  cell_id <- as.integer(cell_id)
  if (!length(cell_id) || anyNA(cell_id) || anyDuplicated(cell_id)) {
    .fire_stop('cell_id must contain unique non-missing identifiers.')
  }
  if (length(last_fire_year) == 1L) last_fire_year <- rep(last_fire_year, length(cell_id))
  if (length(last_severity) == 1L) last_severity <- rep(last_severity, length(cell_id))
  if (length(last_fire_year) != length(cell_id) ||
      length(last_severity) != length(cell_id)) {
    .fire_stop('state vectors must have length 1 or length(cell_id).')
  }
  x <- list(
    cells = data.frame(
      cell_id = cell_id,
      last_fire_year = as.integer(last_fire_year),
      last_severity = as.numeric(last_severity)
    ),
    extra = extra
  )
  class(x) <- 'landisR_fire_state'
  validate_fire_state(x)
  x
}

#' Validate external fire-model state
#' @param x Object to validate.
#' @export
validate_fire_state <- function(x) {
  if (!inherits(x, 'landisR_fire_state')) {
    .fire_stop('x must inherit from landisR_fire_state.')
  }
  if (!is.data.frame(x$cells) ||
      !all(c('cell_id','last_fire_year','last_severity') %in% names(x$cells))) {
    .fire_stop('invalid fire state cell table.')
  }
  if (anyDuplicated(x$cells$cell_id)) .fire_stop('state cell_id must be unique.')
  invisible(TRUE)
}

#' @export
print.landisR_fire_state <- function(x, ...) {
  cat('<landisR_fire_state>\n')
  cat('  cells : ', nrow(x$cells), '\n', sep = '')
  cat('  burned previously : ', sum(!is.na(x$cells$last_fire_year)), '\n', sep = '')
  invisible(x)
}
