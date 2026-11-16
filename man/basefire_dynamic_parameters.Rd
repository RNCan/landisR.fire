\name{basefire_parameter_schedule}
\alias{basefire_parameter_schedule}
\alias{basefire_update_parameters}
\alias{basefire_step_dynamic}
\title{Schedule and apply dynamic Base Fire region parameters}
\description{Change numeric Base Fire region parameters at specified simulation years, while preserving the current fire state, dynamic region-map configuration and random-number generator.}
\usage{
basefire_parameter_schedule(changes)

basefire_update_parameters(prepared, schedule, current_time)

basefire_step_dynamic(simulation, prepared, parameter_schedule, ..., year = simulation$year)
}
\arguments{
  \item{changes}{Data frame with \code{year}, exactly one of \code{region_name} or \code{map_code}, and one or more parameter columns. An \code{NA} in a parameter column means unchanged.}
  \item{prepared}{Current profile returned by \code{basefire_prepare()} or by a previous \code{basefire_step()}.}
  \item{schedule}{Validated schedule, or input data frame for \code{basefire_parameter_schedule()}.}
  \item{current_time}{Current landisR simulation year, as a nonnegative integer.}
  \item{simulation}{Current \code{landisR} simulation.}
  \item{year}{Current or next succession year.}
  \item{parameter_schedule}{Dynamic parameter schedule.}
  \item{...}{Arguments passed unchanged to \code{basefire_step()}, including output writing arguments.}
}
\details{
Updates are effective at the first simulated fire step for which \code{current_time >= year}. Their values persist until another scheduled update. Supported columns are \code{mean_size}, \code{min_size}, \code{max_size}, \code{ignition_probability} and \code{fire_spread_age} (the Original/Base Fire \eqn{K} parameter). Region names and map codes identify regions; this interface does not change spatial region membership or the map-code identities.

Each schedule may include several regions per year. Duplicate region-year keys and unknown regions are errors. The same schedule must be reused throughout one prepared run, and simulation time must not move backwards. When there are no parameter changes, the existing prepared caches are reused. When values change, the entire Base Fire preparation is rebuilt with the original \code{fire_state} and context, preventing stale per-cell caches. An audit table is stored in \code{prepared$context$parameter_change_log}. No RNG draws are introduced by parameter updating.

Dynamic spatial membership from \code{DynamicFireRegionTable} continues to be handled separately by the existing \code{dynamic_region_map_reader}. Parameter updates are applied before the fire step and therefore before any spatial map replacement scheduled during that step.
}
\value{
\code{basefire_parameter_schedule()} returns an immutable schedule data frame. \code{basefire_update_parameters()} returns an updated prepared profile. \code{basefire_step_dynamic()} returns the same list as \code{basefire_step()}, including its \code{prepared_profile} for the next step.
}
\examples{
schedule <- basefire_parameter_schedule(data.frame(
  year = c(10L, 20L),
  region_name = c("North", "North"),
  ignition_probability = c(0.02, 0.04),
  fire_spread_age = c(60L, 45L)
))
print(schedule)
}
