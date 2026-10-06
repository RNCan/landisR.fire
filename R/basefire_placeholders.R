# Internal placeholders defining the modules to be implemented next.
# They are deliberately not exported until their behavior is validated
# against LANDIS-II Base Fire.

.basefire_not_implemented <- function(component) {
  function(...) {
    .fire_stop('Base Fire ', component, ' is not implemented yet.')
  }
}

.basefire_ignition <- .basefire_not_implemented('ignition')
.basefire_size <- .basefire_not_implemented('size')
.basefire_spread <- .basefire_not_implemented('spread')
.basefire_severity <- .basefire_not_implemented('severity')
.basefire_mortality <- .basefire_not_implemented('mortality')
