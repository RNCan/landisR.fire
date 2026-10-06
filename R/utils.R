.fire_stop <- function(...) stop(..., call. = FALSE)

.assert_scalar_character <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .fire_stop(name, ' must be one non-empty character value.')
  }
  invisible(TRUE)
}

.assert_function <- function(x, name) {
  if (!is.function(x)) .fire_stop(name, ' must be a function.')
  invisible(TRUE)
}

.assert_named_list <- function(x, name) {
  if (!is.list(x)) .fire_stop(name, ' must be a list.')
  invisible(TRUE)
}
