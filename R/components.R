.wrap_component <- function(fun, class_name, prepare = NULL, metadata = list()) {
  .assert_function(fun, 'fun')
  .assert_scalar_character(class_name, 'class_name')

  if (!is.null(prepare)) {
    .assert_function(prepare, 'prepare')
  }
  .assert_named_list(metadata, 'metadata')

  attr(fun, 'fire_component_class') <- class_name
  attr(fun, 'fire_prepare') <- prepare
  attr(fun, 'fire_component_metadata') <- metadata
  fun
}

.as_fire_component <- function(fun, class_name) {
  .assert_function(fun, class_name)

  existing <- attr(fun, 'fire_component_class', exact = TRUE)

  if (is.null(existing)) {
    return(.wrap_component(fun, class_name))
  }

  if (!identical(existing, class_name)) {
    .fire_stop(
      class_name,
      ' component has incompatible fire_component_class: ',
      existing,
      '.'
    )
  }

  prepare <- attr(fun, 'fire_prepare', exact = TRUE)
  if (!is.null(prepare)) {
    .assert_function(prepare, paste0(class_name, ' prepare'))
  }

  fun
}

#' Wrap a custom ignition function
#'
#' @param fun User evaluation function.
#' @param prepare Optional preparation function called once by fire_prepare().
#' @param metadata Optional component metadata.
#' @export
fire_ignition_custom <- function(fun, prepare = NULL, metadata = list()) {
  .wrap_component(fun, 'ignition', prepare, metadata)
}

#' Wrap a custom fire-size function
#'
#' @param fun User evaluation function.
#' @param prepare Optional preparation function called once by fire_prepare().
#' @param metadata Optional component metadata.
#' @export
fire_size_custom <- function(fun, prepare = NULL, metadata = list()) {
  .wrap_component(fun, 'size', prepare, metadata)
}

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
    if (length(n) != 1L || is.na(n) || n < 1L) {
      .fire_stop('n must be >= 1.')
    }
    sample(values, size = n, replace = replace, prob = probabilities)
  }

  attr(fun, 'fire_component_class') <- 'size'
  attr(fun, 'fire_prepare') <- NULL
  attr(fun, 'fire_component_metadata') <- list(distribution = 'empirical')
  attr(fun, 'distribution') <- 'empirical'
  attr(fun, 'values') <- values
  attr(fun, 'probabilities') <- probabilities
  attr(fun, 'cdf') <- cumsum(probabilities)
  fun
}

#' Wrap a custom spread function
#'
#' @param fun User evaluation function.
#' @param prepare Optional preparation function called once by fire_prepare().
#' @param metadata Optional component metadata.
#' @export
fire_spread_custom <- function(fun, prepare = NULL, metadata = list()) {
  .wrap_component(fun, 'spread', prepare, metadata)
}

#' Wrap a custom conditional burn-probability function
#'
#' This component answers a different question from spread: spread determines
#' which cells are exposed/reachable, whereas burn_probability determines the
#' probability that an exposed candidate cell actually burns.
#'
#' @param fun User evaluation function.
#' @param prepare Optional preparation function called once by fire_prepare().
#' @param metadata Optional component metadata.
#' @export
fire_burn_probability_custom <- function(fun, prepare = NULL, metadata = list()) {
  .wrap_component(fun, 'burn_probability', prepare, metadata)
}

#' Wrap a custom severity function
#'
#' @param fun User evaluation function.
#' @param prepare Optional preparation function called once by fire_prepare().
#' @param metadata Optional component metadata.
#' @export
fire_severity_custom <- function(fun, prepare = NULL, metadata = list()) {
  .wrap_component(fun, 'severity', prepare, metadata)
}

#' Wrap a custom mortality function
#'
#' @param fun User evaluation function.
#' @param prepare Optional preparation function called once by fire_prepare().
#' @param metadata Optional component metadata.
#' @export
fire_mortality_custom <- function(fun, prepare = NULL, metadata = list()) {
  .wrap_component(fun, 'mortality', prepare, metadata)
}
