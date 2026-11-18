\name{landisR_fire_api}
\alias{landisR_fire_api}
\alias{fire_burn_probability_custom}
\alias{fire_event}
\alias{fire_grid_topology}
\alias{fire_ignition_custom}
\alias{fire_model}
\alias{fire_mortality_custom}
\alias{fire_prepare}
\alias{fire_severity_custom}
\alias{fire_size_custom}
\alias{fire_size_empirical}
\alias{fire_spread_custom}
\alias{fire_state}
\alias{validate_fire_event}
\alias{validate_fire_model}
\alias{validate_fire_prepared_model}
\alias{validate_fire_state}
\alias{validate_fire_topology}
\title{Construct, prepare and validate modular fire models}
\description{
Build a custom fire model from six replaceable components. Prepare a model
against a specific grid to cache calculations, and use the event and state
constructors to represent individual fires and persistent fire history.
}
\usage{
fire_model(ignition, size, spread, burn_probability, severity, mortality,
  name = "custom", metadata = list())
fire_prepare(model, landscape, fire_state = NULL, context = list(), topology = NULL)

fire_grid_topology(landscape, cell_id = "cell_id", row = "row",
  col = "col", active = "active")
fire_state(cell_id, last_fire_year = NA_integer_, last_severity = NA_real_,
  extra = list())
fire_event(year, ignition_cell, cells, severity, target_area = NA_real_,
  actual_area = NA_real_, fire_region = NA_character_, source_cell = NULL,
  metadata = list())

fire_ignition_custom(fun, prepare = NULL, metadata = list())
fire_size_custom(fun, prepare = NULL, metadata = list())
fire_spread_custom(fun, prepare = NULL, metadata = list())
fire_burn_probability_custom(fun, prepare = NULL, metadata = list())
fire_severity_custom(fun, prepare = NULL, metadata = list())
fire_mortality_custom(fun, prepare = NULL, metadata = list())
fire_size_empirical(values, probabilities, replace = TRUE)

validate_fire_model(x)
validate_fire_prepared_model(x)
validate_fire_topology(x)
validate_fire_state(x)
validate_fire_event(x)
}
\arguments{
  \item{ignition, size, spread, burn_probability, mortality}{Functions or
  corresponding wrapped components. They determine ignitions, target area,
  cells exposed to fire, conditional probability of burning, and cohort
  mortality, respectively.}
  \item{severity}{In \code{fire_model()}, a function or wrapped component
  that assigns severity. In \code{fire_event()}, a finite scalar, vector
  aligned to \code{cells}, or table with matching \code{cell_id} and
  \code{severity} columns.}
  \item{name}{Display name of the fire model.}
  \item{metadata}{Named list of model, component, or event metadata.}
  \item{model}{A \code{landisR_fire_model}.}
  \item{landscape}{For topology, a data frame with IDs, rows and columns.
  For model preparation, the landscape object supplied to the callbacks.}
  \item{fire_state}{Optional fire-history object supplied to preparers.}
  \item{context}{Named list of additional context supplied to preparers.}
  \item{topology}{Optional prebuilt topology. If omitted, the preparer
  builds one when the default grid columns are present.}
  \item{cell_id}{For topology, the name of the ID column. For fire state,
  a vector of unique cell IDs.}
  \item{row, col, active}{Column names for grid coordinates and the
  optional logical active mask. All rows are active if the mask is absent.}
  \item{last_fire_year, last_severity}{Scalar or per-cell fire-history
  vectors. Use missing values for unknown history.}
  \item{extra}{Named list of additional state fields.}
  \item{year}{Simulation year of the event.}
  \item{ignition_cell}{One cell ID included in \code{cells}.}
  \item{cells}{Unique burned cell IDs in spread order.}
  \item{target_area, actual_area}{Optional target and realized fire
  areas in a consistent unit.}
  \item{fire_region}{Optional identifier of the ignition region.}
  \item{source_cell}{Optional predecessor ID per burned cell.}
  \item{fun}{Evaluation function to tag as the given component type.
  Its arguments and return value are determined by the calling profile.}
  \item{prepare}{Optional function accepting the named arguments
  \code{landscape}, \code{fire_state}, \code{topology}, \code{context}
  and \code{model}, either explicitly or through \code{...}.
  Its return value is the cached component data.}
  \item{values}{Finite positive numeric candidate fire sizes, typically
  in hectares.}
  \item{probabilities}{Nonnegative weights, one per candidate size.
  At least one must be positive. Weights are normalized to sum to one.}
  \item{replace}{Sample the candidate sizes with replacement.}
  \item{x}{Object of the type named by the validator.}
}
\details{
\code{fire_grid_topology()} stores internal neighbor indices for two
directional rings. Index zero denotes an absent or inactive neighbor.
\code{fire_prepare()} runs each preparation callback once. The resulting
caches are available in the prepared model's \code{prepared} list.
For a strict Base Fire run, see \code{\link{basefire_prepare}}.
}
\value{
\code{fire_model()}, \code{fire_prepare()}, \code{fire_grid_topology()},
\code{fire_state()} and \code{fire_event()} return, respectively,
\code{landisR_fire_model}, \code{landisR_fire_prepared_model},
\code{landisR_fire_topology}, \code{landisR_fire_state} and
\code{landisR_fire_event} objects. An event contains a cell-level severity
table and a spread table with predecessors. Custom component constructors
return the supplied function with its component attributes attached.
\code{fire_size_empirical()} returns a callable size component that samples
\code{n} values (default one). Validators return \code{TRUE} invisibly on success and
signal an error for invalid objects.
}
\seealso{\code{\link{basefire_read}}, \code{\link{basefire_prepare}},
  \code{\link{basefire_step}}}
\examples{
grid <- data.frame(cell_id = 1:3, row = 1L, col = 1:3)
topology <- fire_grid_topology(grid)
state <- fire_state(grid$cell_id)
event <- fire_event(10L, ignition_cell = 2L, cells = c(2L, 3L),
  severity = c(3, 2), source_cell = c(NA_integer_, 2L))
validate_fire_topology(topology)
validate_fire_state(state)
validate_fire_event(event)

# Assemble a toy profile and cache its neighbor matrix.
ignition <- fire_ignition_custom(function(...) 2L)
size <- fire_size_custom(function(...) 2)
spread <- fire_spread_custom(
  function(...) c(2L, 3L),
  prepare = function(topology, ...) topology$neighbors
)
burn <- fire_burn_probability_custom(function(...) 1)
severity <- fire_severity_custom(function(...) 3)
mortality <- fire_mortality_custom(function(cohorts, ...) cohorts$age > 10L)
model <- fire_model(ignition, size, spread, burn, severity, mortality)
validate_fire_model(model)
prepared <- fire_prepare(model, grid, state, topology = topology)
validate_fire_prepared_model(prepared)

draw_size <- fire_size_empirical(c(10, 100, 500), c(6, 3, 1))
draw_size(n = 3L)
}
\concept{fire}
\concept{modular fire models}
