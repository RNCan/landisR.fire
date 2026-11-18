.fire_elapsed <- function() unname(proc.time()[["elapsed"]])

.fire_timing_enabled <- function() is.environment(getOption("landisR.fire.timing"))

.fire_add_timing <- function(level, year, phases, seconds) {
  sink <- getOption("landisR.fire.timing")
  if (!is.environment(sink)) return(invisible(NULL))
  row <- data.frame(annee = as.integer(year), niveau = level, phase = phases,
                    secondes = as.numeric(seconds))
  sink$phases[[length(sink$phases) + 1L]] <- row
  invisible(NULL)
}

.fire_timer <- function(level, year) {
  if (!.fire_timing_enabled()) return(function(phase) invisible(NULL))
  force(level)
  force(year)
  last <- .fire_elapsed()
  function(phase) {
    now <- .fire_elapsed()
    .fire_add_timing(level, year, phase, now - last)
    last <<- .fire_elapsed()
    invisible(NULL)
  }
}
