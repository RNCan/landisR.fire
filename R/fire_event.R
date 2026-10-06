#' Construct a standard fire event
#'
#' The event object is the exchange format between a fire generator and the
#' forest-impact/application layer. Severity is explicitly cell based.
#'
#' @param year Event year.
#' @param ignition_cell Ignition cell identifier.
#' @param cells Burned cell identifiers in fire-spread order.
#' @param severity Either one scalar severity, one value per cell, or a
#'   data.frame with cell_id and severity columns.
#' @param target_area Optional target fire area.
#' @param actual_area Optional realized fire area.
#' @param fire_region Optional ignition-region identifier.
#' @param source_cell Optional source/predecessor cell for each burned cell.
#' @param metadata Optional event metadata.
#' @export
fire_event <- function(year, ignition_cell, cells, severity,
                       target_area = NA_real_, actual_area = NA_real_,
                       fire_region = NA_character_, source_cell = NULL,
                       metadata = list()) {
  cells <- as.integer(cells)
  if (!length(cells) || anyNA(cells) || anyDuplicated(cells)) {
    .fire_stop('cells must contain unique non-missing cell identifiers.')
  }
  ignition_cell <- as.integer(ignition_cell)
  if (length(ignition_cell) != 1L || is.na(ignition_cell) ||
      !ignition_cell %in% cells) {
    .fire_stop('ignition_cell must be one member of cells.')
  }

  if (is.data.frame(severity)) {
    if (!all(c('cell_id','severity') %in% names(severity))) {
      .fire_stop('severity data.frame must contain cell_id and severity.')
    }
    severity_df <- severity[, c('cell_id','severity'), drop = FALSE]
    severity_df$cell_id <- as.integer(severity_df$cell_id)
    if (anyDuplicated(severity_df$cell_id) ||
        !setequal(severity_df$cell_id, cells)) {
      .fire_stop('severity cell_id values must match cells exactly.')
    }
    severity_df <- severity_df[match(cells, severity_df$cell_id), , drop = FALSE]
  } else {
    severity <- as.numeric(severity)
    if (length(severity) == 1L) severity <- rep(severity, length(cells))
    if (length(severity) != length(cells) || any(!is.finite(severity))) {
      .fire_stop('severity must be finite and scalar or one value per cell.')
    }
    severity_df <- data.frame(cell_id = cells, severity = severity)
  }

  if (is.null(source_cell)) source_cell <- rep(NA_integer_, length(cells))
  source_cell <- as.integer(source_cell)
  if (length(source_cell) != length(cells)) {
    .fire_stop('source_cell must be NULL or one value per burned cell.')
  }

  x <- list(
    year = as.integer(year),
    ignition_cell = ignition_cell,
    cells = cells,
    severity = severity_df,
    target_area = as.numeric(target_area),
    actual_area = as.numeric(actual_area),
    fire_region = as.character(fire_region),
    spread = data.frame(
      order = seq_along(cells),
      cell_id = cells,
      source_cell = source_cell
    ),
    metadata = metadata
  )
  class(x) <- 'landisR_fire_event'
  validate_fire_event(x)
  x
}

#' Validate a standard fire event
#' @param x Object to validate.
#' @export
validate_fire_event <- function(x) {
  if (!inherits(x, 'landisR_fire_event')) {
    .fire_stop('x must inherit from landisR_fire_event.')
  }
  if (!length(x$cells) || anyDuplicated(x$cells)) {
    .fire_stop('fire event cells must be unique.')
  }
  if (!x$ignition_cell %in% x$cells) {
    .fire_stop('ignition cell must belong to event cells.')
  }
  if (!setequal(x$severity$cell_id, x$cells)) {
    .fire_stop('severity map must match event cells.')
  }
  invisible(TRUE)
}

#' @export
print.landisR_fire_event <- function(x, ...) {
  cat('<landisR_fire_event>\n')
  cat('  year          : ', x$year, '\n', sep = '')
  cat('  ignition cell : ', x$ignition_cell, '\n', sep = '')
  cat('  burned cells  : ', length(x$cells), '\n', sep = '')
  cat('  severity range: ', paste(range(x$severity$severity), collapse = ' - '), '\n', sep = '')
  invisible(x)
}
