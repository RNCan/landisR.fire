\name{landisR.fire-package}
\alias{landisR.fire-package}
\alias{landisR.fire}
\title{landisR.fire: Fire Disturbance for LANDIS-R}
\description{
Run the Base Fire compatibility profile on a shared LANDIS-R forest or
assemble a modular fire model for a custom profile.
}
\details{
Read parameters with \code{\link{basefire_read}} or construct them in R
with \code{\link{basefire_build}}. Register \code{\link{basefire_extension}}
in a coupled simulation. For a manual loop, use
\code{\link{basefire_prepare}} and \code{\link{basefire_step}}.
}
\seealso{
\code{\link{basefire_mortality}}, \code{\link{basefire_write_outputs}},
\code{\link{basefire_parameter_schedule}}, \code{\link{fire_model}}
}
\docType{package}
