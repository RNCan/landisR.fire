.fire_component_names <- c(
  'ignition',
  'size',
  'spread',
  'burn_probability',
  'severity',
  'mortality'
)

#' Construct a modular fire model
#'
#' A fire model is an ordered collection of replaceable process components.
#' Base Fire will later be represented by one strict combination of these
#' components rather than by a monolithic implementation.
#'
#' @param ignition Function that proposes ignition events.
#' @param size Function that samples or computes target fire size.
#' @param spread Function that determines which cells are spatially exposed.
#' @param burn_probability Function that determines the conditional probability
#'   that an exposed candidate cell burns.
#' @param severity Function that assigns cell-level fire severity.
#' @param mortality Function that converts fire exposure/severity to forest impact.
#' @param name Human-readable model name.
#' @param metadata Optional model metadata.
#' @export
fire_model <- function(ignition, size, spread, burn_probability,
                       severity, mortality, name = 'custom', metadata = list()) {
  components <- list(
    ignition = .as_fire_component(ignition, 'ignition'),
    size = .as_fire_component(size, 'size'),
    spread = .as_fire_component(spread, 'spread'),
    burn_probability = .as_fire_component(
      burn_probability,
      'burn_probability'
    ),
    severity = .as_fire_component(severity, 'severity'),
    mortality = .as_fire_component(mortality, 'mortality')
  )

  .assert_scalar_character(name, 'name')
  .assert_named_list(metadata, 'metadata')

  x <- c(
    list(name = name),
    components,
    list(metadata = metadata)
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

  for (nm in .fire_component_names) {
    .assert_function(x[[nm]], nm)

    component_class <- attr(
      x[[nm]],
      'fire_component_class',
      exact = TRUE
    )

    if (!identical(component_class, nm)) {
      .fire_stop(
        nm,
        ' component must have fire_component_class="',
        nm,
        '".'
      )
    }
  }

  invisible(TRUE)
}

#' @export
print.landisR_fire_model <- function(x, ...) {
  cat('<landisR_fire_model>\n')
  cat('  name       : ', x$name, '\n', sep = '')
  cat(
    '  components : ',
    paste(.fire_component_names, collapse = ', '),
    '\n',
    sep = ''
  )
  invisible(x)
}
