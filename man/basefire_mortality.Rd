\name{basefire_mortality}
\alias{basefire_mortality}
\title{Report cohort and biomass mortality from Base Fire}
\usage{basefire_mortality(result, extension = "fire", by_species = FALSE)}
\arguments{
\item{result}{A coupled run from \code{landisR::landis_sim()} or a single
  step from \code{basefire_step()}.}
\item{extension}{Name of the Base Fire disturbance in a coupled run.}
\item{by_species}{If \code{TRUE}, return one row per species killed by a
  logged fire.}
}
\value{A data frame with one row per logged fire, the existing Base Fire
log fields including \code{cohorts_killed}, and \code{MgBiomassKilled},
the total living biomass killed in megagrams (metric tonnes).
With \code{by_species = TRUE}, columns are \code{time}, \code{event_id},
\code{species}, \code{cohorts_killed} and \code{MgBiomassKilled}. If the
historical species detail is unavailable, the event total is reported in
one row with \code{species = NA}.}
\description{Reports the biomass of cohorts immediately before removal. With
Biomass Succession, cohort biomass in grams per square metre is converted
using the physical project cell area. A fire with no killed cohorts has
zero killed biomass. A run without Biomass Succession, or a historical
coupled event with deaths whose biomass was not retained, returns
\code{NA} for unavailable biomass.}
\details{The reported quantity is live cohort biomass killed. It does not
describe biomass consumed or emissions. The event identifier is local to
a simulation year. For coupled runs, use the year and event identifier
together to identify a fire.}
\seealso{\code{basefire_step()}, \code{basefire_extension()}, \code{basefire_write_outputs()}}
\examples{
\dontrun{
mortality <- basefire_mortality(run, extension = "fire")
write.csv(mortality, "fire_mortality.csv", row.names = FALSE)
mortality_species <- basefire_mortality(run, by_species = TRUE)
}
}
