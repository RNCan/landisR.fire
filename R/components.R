.wrap_component <- function(fun, class_name) {
  .assert_function(fun, 'fun')
  attr(fun, 'fire_component_class') <- class_name
  fun
}

#' Wrap a custom ignition function
#' @param fun User function.
#' @export
fire_ignition_custom <- function(fun) .wrap_component(fun, 'ignition')

#' Wrap a custom fire-size function
#' @param fun User function.
#' @export
fire_size_custom <- function(fun) .wrap_component(fun, 'size')

#' Empirical/custom-frequency fire-size distribution
#'
#' @param values Possible fire sizes, normally in hectares or another unit
#'   chosen consistently by the caller.
#' @param probabilities Sampling probabilities associated with values.
#' @param replace Passed to sample().
#' @export
fire_size_empirical <- function(values, probabilities, replace = TRUE) {
  if (!is.numeric(values) || !length(values) || any(!is.finite(values)) ||
      any(values <= 0)) {
    .fire_stop('values must contain finite positive numeric fire sizes.')
  }
  if (!is.numeric(probabilities) || length(probabilities) != length(values) ||
      any(!is.finite(probabilities)) || any(probabilities < 0) ||
      sum(probabilities) <= 0) {
    .fire_stop('probabilities must be non-negative and match values.')
  }
  probabilities <- probabilities / sum(probabilities)

  fun <- function(n = 1L, ...) {
    n <- as.integer(n)
    if (length(n) != 1L || is.na(n) || n < 1L) .fire_stop('n must be >= 1.')
    sample(values, size = n, replace = replace, prob = probabilities)
  }
  attr(fun, 'fire_component_class') <- 'size'
  attr(fun, 'distribution') <- 'empirical'
  attr(fun, 'values') <- values
  attr(fun, 'probabilities') <- probabilities
  fun
}

#' Wrap a custom spread function
#' @param fun User function.
#' @export
fire_spread_custom <- function(fun) .wrap_component(fun, 'spread')

#' Wrap a custom severity function
#' @param fun User function.
#' @export
fire_severity_custom <- function(fun) .wrap_component(fun, 'severity')

#' Wrap a custom mortality function
#' @param fun User function.
#' @export
fire_mortality_custom <- function(fun) .wrap_component(fun, 'mortality')
