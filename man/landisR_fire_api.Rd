\name{landisR_fire_api}
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
\title{Core modular fire-model API}
\description{
Core constructors, validators, preparation utilities, topology helpers,
event and state objects, and customizable process components used by
\pkg{landisR.fire}. These functions provide the generic modular fire-model
interface on which specific fire profiles, including the strict
Base Fire compatibility profile, are built.
}
\details{
The generic fire API separates fire processes into replaceable components
for ignition, target fire size, spread, conditional burn probability,
severity, and mortality. A fire model can be prepared once against a
landscape, external fire state, topology, and model context, allowing
component-specific caches to be reused during simulation.

\code{fire_model()} combines process components into a fire-model object.
\code{fire_prepare()} prepares a model for a specific landscape and state.
\code{fire_state()} stores persistent cell-level fire history.
\code{fire_event()} stores one realized fire event and its cell-level
severity information.
\code{fire_grid_topology()} constructs active-cell topology from a
landscape grid.

The \code{fire_*_custom()} constructors wrap user-supplied functions as
validated fire components. \code{fire_size_empirical()} provides the
built-in empirical fire-size component.

The \code{validate_*()} functions validate the corresponding fire objects
and return invisibly when validation succeeds.
}
\value{
Constructors return objects representing fire models, prepared models,
components, topology, state, or events, depending on the function called.
Validation functions return invisibly on success and signal an error when
the supplied object is invalid.
}
\seealso{
\code{fire_model}, \code{fire_prepare}, \code{fire_state},
\code{fire_event}, \code{basefire_read}, \code{basefire_prepare},
\code{basefire_step}
}
\concept{fire}
\concept{disturbance}
\concept{modular fire models}
