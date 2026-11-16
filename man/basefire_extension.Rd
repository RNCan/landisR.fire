\name{basefire_extension}
\alias{basefire_extension}
\title{Run Base Fire before succession in landis_sim}
\usage{
basefire_extension(parameters, fire_region_map_code, cell_area = NULL,
  fire_state = NULL, dynamic_region_map_reader = NULL,
  time_of_last_wind = NULL, parameter_schedule = NULL, seed = NULL,
  stochastic_serotiny = TRUE, rng_mode = c("independent", "landis_v7"))
}
\arguments{
  \item{parameters, fire_region_map_code, fire_state, dynamic_region_map_reader,
  time_of_last_wind}{Base Fire settings passed to \code{basefire_prepare()}.}
  \item{cell_area}{Cell area in hectares. With \code{NULL}, derive the area
  from the project cell length in metres.}
  \item{parameter_schedule}{Optional schedule for numeric fire-region parameters.}
  \item{seed}{A nonnegative integer for the fire engine's R random-number
  stream. With \code{NULL}, initialize from the simulation's R RNG state.
  Use \code{NULL} with \code{rng_mode = "landis_v7"}.}
  \item{stochastic_serotiny}{Control stochastic serotinous recruitment.}
  \item{rng_mode}{With \code{"independent"} (default), initialize a separate
  R stream for the fire engine. With \code{"landis_v7"}, fire ignition, size
  and spread use the simulation's shared LANDIS-II v7 MT19937 stream, in
  the configured order of disturbances. This mode requires the native
  fire backend and a simulation initialized with
  \code{recruitment_rng = "landis_v7"}.}
}
\description{
Creates an extension for \code{landisR::landis_sim()}. It runs in declaration
order with other disturbances before succession. Harvest reads the sparse
fire-date history on the shared forest state. Severity-positive cells enter
the history even when the event causes no cohort mortality.
}
\value{
A disturbance descriptor. At each step it returns a sparse \code{burned}
table (\code{cell_id} and \code{severity}), a fire event \code{log} and a
\code{summary}. The prepared profile and fire-engine R RNG state are kept
with the coupled run for checkpointing. The shared MT19937 stream is kept
with the simulation, including after fire timesteps with no ignitions.
}
