\name{basefire}
\alias{basefire_read}
\alias{basefire_prepare}
\alias{basefire_step}
\title{Run the LANDIS-II Base Fire compatibility profile with landisR}
\description{
Read LANDIS-II Base Fire parameters, prepare the strict Base Fire
compatibility profile against a live \pkg{landisR} simulation, and execute
one fire timestep while applying mortality and post-fire recruitment through
the public \pkg{landisR} state/update/fire interface.
}
\usage{
basefire_read(path, known_species = NULL, base_dir = dirname(path))

basefire_prepare(
  simulation,
  parameters,
  fire_region_map_code,
  cell_area,
  fire_state = NULL,
  dynamic_region_map_reader = NULL,
  time_of_last_wind = NULL
)

basefire_step(
  simulation,
  prepared,
  rng = list(),
  stochastic_serotiny = TRUE,
  metadata = list(), progress = FALSE, verbose = FALSE,
  output_dir = NULL, output_formats = c("csv", "rds"),
  raster_template = NULL, cell_length = NULL, overwrite_outputs = FALSE
)
}
\arguments{
  \item{output_dir}{Optional parent directory for per-timestep fire outputs. NULL disables writing.}
  \item{output_formats}{Subset of \code{c("csv", "rds", "tif")}; default CSV and RDS.}
  \item{raster_template}{Optional georeferenced one-layer terra raster (or path) used for GeoTIFF grid geometry.}
  \item{cell_length}{Optional relative grid cell size if no georeferenced raster is supplied.}
  \item{overwrite_outputs}{Allow replacing an existing annual output folder; default FALSE.}
  \item{progress}{Logical. Show progress information in the console; default \code{FALSE}.}
  \item{verbose}{Logical. Show detailed console messages; default \code{FALSE}.}
  \item{path}{Path to a LANDIS-II Base Fire parameter file.}
  \item{known_species}{Optional character vector of species codes allowed
  in the Base Fire species table.}
  \item{base_dir}{Base directory used to resolve relative paths referenced
  by the Base Fire parameter file.}
  \item{simulation}{A live \code{landisR_simulation}.}
  \item{parameters}{A Base Fire parameter object, normally returned by
  \code{basefire_read()}.}
  \item{fire_region_map_code}{Fire-region map codes aligned to the landisR
  landscape. May be a scalar, a named vector keyed by cell identifier, or a
  data frame with \code{cell_id} and \code{fire_region_map_code}.}
  \item{cell_area}{Positive cell area expressed in the same area units as
  the Base Fire mean, minimum, and maximum fire-size parameters.}
  \item{fire_state}{Optional persistent \code{landisR_fire_state}. If
  omitted, initial time since fire is inferred from maximum cohort age.}
  \item{dynamic_region_map_reader}{Optional function used to supply
  scheduled dynamic Base Fire region maps.}
  \item{time_of_last_wind}{Optional per-cell time of last wind disturbance
  used by Base Fire severity calculations.}
  \item{prepared}{Prepared Base Fire profile returned by
  \code{basefire_prepare()} or by the previous \code{basefire_step()}.}
  \item{rng}{Named list of optional Base Fire random-number controls used
  primarily for deterministic testing or specialized simulation control.}
  \item{stochastic_serotiny}{Logical value passed to
  \code{landisR::landis_fire()} to control stochastic post-fire
  serotinous recruitment.}
  \item{metadata}{Named list appended to the landisR fire-event metadata.}
}
\details{
\code{basefire_read()} parses the strict LANDIS-II Base Fire parameter
format and preserves the compatibility rules implemented by
\pkg{landisR.fire}.

\code{basefire_prepare()} binds those parameters to a live
\pkg{landisR} simulation. The current coupled adapter requires the Base
Fire timestep to equal the landisR succession timestep.

\code{basefire_step()} reads the current cohort state with
\code{landisR::landis_state()}, runs the strict Base Fire event chain,
converts post-fire survivors and completely emptied cells to a
\code{landisR::landis_update()}, and applies that update through
\code{landisR::landis_fire()}. This ordering allows landisR to capture
pre-fire serotinous seed sources before mortality is applied and then perform
post-fire recruitment.

The function does not advance succession with \code{landisR::landis_step()}.
The caller controls succession timing explicitly. For repeated fire
timesteps, carry both \code{result$simulation} and
\code{result$prepared_profile} forward.
}
\value{
\code{basefire_read()} returns validated strict Base Fire parameters.

\code{basefire_prepare()} returns a prepared strict Base Fire profile.

\code{basefire_step()} returns a list with:
\itemize{
  \item \code{simulation}: the updated \code{landisR_simulation};
  \item \code{fire_result}: the complete strict Base Fire timestep result;
  \item \code{prepared_profile}: the profile carrying updated fire history;
  \item \code{update}: the applied \code{landisR_update}, or \code{NULL}
  when no cell burned;
  \item \code{burned_cells}: integer identifiers of burned cells.
}
}
\seealso{
\code{fire_state}, \code{fire_event}, \code{fire_model}
}
\concept{fire}
\concept{Base Fire}
\concept{landisR coupling}
