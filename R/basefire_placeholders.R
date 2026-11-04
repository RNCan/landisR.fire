# Internal placeholders defining the Base Fire modules to be implemented next.
# They are deliberately not exported until their behavior is validated
# against LANDIS-II Base Fire.
#
# Ignition is implemented in basefire_ignition.R.
# Target fire size is implemented in basefire_size.R.
# Wind-driven candidate spread and the TSF burn-probability gate are
# implemented in basefire_spread.R.

.basefire_not_implemented <- function(component) {
  function(...) {
    .fire_stop('Base Fire ', component, ' is not implemented yet.')
  }
}

.basefire_severity <- .basefire_not_implemented('severity')
.basefire_mortality <- .basefire_not_implemented('mortality')
