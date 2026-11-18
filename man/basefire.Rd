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
basefire_read(path, known_species = NULL, base_dir = dirname(path),
  project = NULL, scenario_directory = dirname(path), raster_template = NULL)

basefire_prepare(
  simulation,
  parameters,
  fire_region_map_code = NULL,
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
  metadata = list(), progress = getOption("landisR.progress", interactive()),
  verbose = FALSE,
  output_dir = NULL, output_formats = c("csv", "rds"),
  raster_template = NULL, cell_length = NULL, overwrite_outputs = FALSE,
  year = simulation$year
)
}
\arguments{
  \item{output_dir}{Optional parent directory for per-timestep fire outputs. NULL disables writing.}
  \item{output_formats}{Subset of \code{c("csv", "rds", "tif")}; default CSV and RDS.}
  \item{raster_template}{For \code{basefire_read()}, an optional path to a
  one-layer landscape raster used to verify the fire maps. Relative paths
  use \code{scenario_directory}. The initial fire map provides the reference
  grid when this argument is omitted. For output, a one-layer terra raster
  or a path defines the GeoTIFF grid geometry.}
  \item{scenario_directory}{Directory used to resolve relative fire raster paths.
  Set this to the directory containing the original LANDIS-II scenario when
  maps are referenced relative to that scenario.}
  \item{cell_length}{Optional relative grid cell size if no georeferenced raster is supplied.}
  \item{overwrite_outputs}{Allow replacing an existing annual output folder; default FALSE.}
  \item{progress}{Show a console bar for the fire step and optional output
  writing. Defaults to display in interactive sessions. For coupled runs,
  use the progress argument of \code{landisR::landis_sim()}.}
  \item{verbose}{Logical. Show detailed console messages; default \code{FALSE}.}
  \item{path}{Path to a LANDIS-II Base Fire parameter file.}
  \item{known_species}{Optional character vector of species codes allowed
  in the Base Fire species table.}
  \item{base_dir}{Directory used to resolve a relative species CSV path.}
  \item{project}{A \code{landisR_project}, typically from \code{landis_read()}.
  For an original LANDIS-II file, the project provides the species codes and
  fire tolerances. A \code{Species_CSV_File} entry explicitly selects the CSV
  specified in the file.}
  \item{year}{Current or next succession year. The next year is used for a
  fire applied before succession.}
  \item{simulation}{A live \code{landisR_simulation}.}
  \item{parameters}{A Base Fire parameter object, normally returned by
  \code{basefire_read()} or \code{basefire_build()}.}
  \item{fire_region_map_code}{Fire-region map codes aligned to the landisR
  landscape. May be a scalar, a named vector keyed by cell identifier, or a
  data frame with \code{cell_id} and \code{fire_region_map_code}.
  When \code{NULL}, use the map stored by \code{basefire_build()}.}
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
The native fused Base Fire backend is selected by default for the strict
model with standard random-number and wind controls. Custom RNG or wind
hooks use the legacy backend when no backend option is set. Explicitly set
\code{options(landisR.fire.backend = "legacy")} to use the legacy engine,
or \code{options(landisR.fire.backend = "native")} to require the native
engine. Explicit native selection rejects unsupported RNG or wind hooks.
The backend option does not change the succession engine.

\code{basefire_read()} imports the original LANDIS-II Base Fire parameter
file with fire tolerances from \code{project}. Inputs with a
\code{Species_CSV_File} entry remain supported.
The returned parameters can be passed directly to \code{basefire_extension()}
or \code{basefire_prepare()}. With \pkg{terra}, the initial raster is read
at preparation and scheduled rasters are read when their year is reached.
A supplied \code{fire_region_map_code} or \code{dynamic_region_map_reader}
takes precedence over the corresponding automatic read.

\code{basefire_prepare()} binds those parameters to a live
\pkg{landisR} simulation. The current coupled adapter requires the Base
Fire timestep to equal the landisR succession timestep.

\code{basefire_step()} reads the shared cohort state, runs the strict Base Fire event chain,
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

\examples{
\dontrun{
# Use an existing canonical biomass project and Base Fire input files.
project <- landisR::landis_from_canonical("inputs/canonical")
sim <- landisR::landis_initialize(project)
parameters <- basefire_read(
  "inputs/base-fire.txt",
  project = project,
  scenario_directory = "inputs",
  raster_template = "ecoregions.tif"
)
area_ha <- project$canonical$cell_length^2 / 10000
prepared <- basefire_prepare(sim, parameters, cell_area = area_ha)
step <- basefire_step(sim, prepared, year = project$canonical$timestep,
  progress = FALSE)
sim <- landisR::landis_step(step$simulation, progress = FALSE)
prepared <- step$prepared_profile
}
}
