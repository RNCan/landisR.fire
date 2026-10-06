#' Construct a modular fire model
#'
#' A fire model is an ordered collection of replaceable process components.
#' Base Fire will later be represented by one strict combination of these
#' components rather than by a monolithic implementation.
#'
#' @param ignition Function that proposes ignition events.
#' @param size Function that samples or computes target fire size.
#' @param spread Function that propagates a fire across cells.
#' @param severity Function that assigns cell-level fire severity.
#' @param mortality Function that converts fire exposure/severity to forest impact.
#' @param name Human-readable model name.
#' @param metadata Optional model metadata.
#' @export
fire_model <- function(ignition, size, spread, severity, mortality,
                       name = 'custom', metadata = list()) {
  .assert_function(ignition, 'ignition')
  .assert_function(size, 'size')
  .assert_function(spread, 'spread')
  .assert_function(severity, 'severity')
  .assert_function(mortality, 'mortality')
  .assert_scalar_character(name, 'name')
  .assert_named_list(metadata, 'metadata')

  x <- list(
    name = name,
    ignition = ignition,
    size = size,
    spread = spread,
    severity = severity,
    mortality = mortality,
    metadata = metadata
  )
  class(x) <- 'landisR_fire_model'
  validate_fire_model(x)
  x
}

#' Validate a modular fire model
#' @param x Object to validate.
#' @export
validate_fire_model <- function(x) {
  if (!inherits(x, 'landisR_fire_model')) {
    .fire_stop('x must inherit from landisR_fire_model.')
  }
  for (nm in c('ignition','size','spread','severity','mortality')) {
    .assert_function(x[[nm]], nm)
  }
  invisible(TRUE)
}

#' @export
print.landisR_fire_model <- function(x, ...) {
  cat('<landisR_fire_model>\n')
  cat('  name       : ', x$name, '\n', sep = '')
  cat('  components : ignition, size, spread, severity, mortality\n')
  invisible(x)
}
