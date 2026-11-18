\name{basefire_build}
\alias{basefire_build}
\title{Build Base Fire parameters and maps from R tables}
\usage{
basefire_build(timestep, species, regions, damage, fire_region_map_code,
  fuel_curves = NULL, wind_curves = NULL, dynamic_region_maps = NULL,
  map_names_template = "fire-severity-{timestep}.tif",
  log_file = "fire-events.csv", summary_log_file = "fire-summary.csv")
}
\arguments{
\item{timestep}{Positive whole number in years, matching the succession step.}
\item{species}{One row per species with \code{species_code} and integer
  \code{fire_tolerance}. Species codes match the project species table.}
\item{regions}{One row per fire region with \code{name}, \code{map_code},
  \code{mean_size}, \code{min_size}, \code{max_size},
  \code{ignition_probability} and \code{fire_spread_age}.
  Sizes are in hectares when the simulation uses the default cell-area
  conversion. Region order controls the compatible engine's traversal.}
\item{damage}{Ordered damage classes with \code{max_age}, a fraction of
  species longevity, and integer \code{severity_tolerance_diff}.}
\item{fire_region_map_code}{Initial map table with integer \code{cell_id} and
  \code{fire_region_map_code}. Supply one row per landscape cell including
  inactive cells. Active cells must identify a defined fire region.}
\item{fuel_curves, wind_curves}{Optional tables with \code{region_name} and
  \code{severity1} through \code{severity5} in severity order.}
\item{dynamic_region_maps}{Optional long table with \code{year},
  \code{cell_id} and \code{fire_region_map_code}. Each scheduled fire year
  supplies a complete map, including inactive cells.}
\item{map_names_template, log_file, summary_log_file}{Output names stored
  in the parameters for compatible workflows.}
}
\description{Constructs parameters for the strict Base Fire model from
data frames created in R or loaded using \code{read.csv()}.}
\value{A Base Fire parameter object accepted by \code{basefire_extension()}
and \code{basefire_prepare()}. The map and scheduled maps stored in the
object are used when the corresponding arguments are omitted.}
\seealso{\code{basefire_read()}, \code{basefire_parameter_schedule()}}
\examples{
regions <- data.frame(name = "north", map_code = 1L, mean_size = 50,
  min_size = 1, max_size = 100, ignition_probability = 0.005,
  fire_spread_age = 30L)
p <- basefire_build(
  timestep = 10L, species = data.frame(species_code = "PICE.MAR",
    fire_tolerance = 3L),
  regions = regions, damage = data.frame(max_age = 1,
    severity_tolerance_diff = 4L),
  fire_region_map_code = data.frame(cell_id = 1:4, fire_region_map_code = 1L)
)
fire <- basefire_extension(p)
}
