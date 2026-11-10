\name{basefire_write_outputs}
\alias{basefire_write_outputs}
\title{Write Base Fire outputs to annual folders}
\description{Writes tabular fire events, summaries and history, optional sparse RDS event objects, and optional GeoTIFF maps without changing the simulation state or its random number generator.}
\usage{
basefire_write_outputs(result, output_dir, formats = c("csv", "rds"),
  raster_template = NULL, cell_length = NULL, overwrite = FALSE)
}
\arguments{
  \item{result}{A \code{basefire_step} output or \code{fire_result}.}
  \item{output_dir}{Root directory for annual outputs.}
  \item{formats}{One or more of \code{"csv"}, \code{"rds"} and \code{"tif"}.}
  \item{raster_template}{Optional one-layer reference \code{terra::SpatRaster} or path with verified extent and CRS.}
  \item{cell_length}{Optional numeric cell size for relative rasters without verified CRS.}
  \item{overwrite}{Whether to replace an existing annual output directory.}
}
\value{Invisibly, a data frame manifest of files created.}
\details{Each timestep is saved to \code{year_XXXX}. GeoTIFF output requires the optional terra package. Without raster_template, the output grid is in relative coordinates only.}
\examples{
\dontrun{
basefire_write_outputs(step_result, "outputs/fire")
}
}
