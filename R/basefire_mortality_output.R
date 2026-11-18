.basefire_killed_density <- function(result) {
  if (!is.null(result$killed_density)) return(result$killed_density)
  log <- result$event_log
  # Native events already retain the killed cohort rows while the step runs.
  vapply(seq_len(nrow(log)), function(i) {
    if (log$cohorts_killed[i] == 0L) return(0)
    id <- log$event_id[i]
    if (id > length(result$raw_events)) return(NA_real_)
    density <- result$raw_events[[id]]$affected_cohorts$biomass
    if (is.null(density)) NA_real_ else sum(density)
  }, numeric(1L))
}

.basefire_killed_species <- function(result) {
  if (!is.null(result$killed_species)) return(result$killed_species)
  lapply(result$event_log$event_id, function(id) {
    if (id > length(result$raw_events)) return(NULL)
    killed <- result$raw_events[[id]]$affected_cohorts
    if (is.null(killed) || !nrow(killed) ||
        is.null(killed$species_code) || anyNA(killed$species_code)) return(NULL)
    n <- nrow(killed)
    density <- if (is.null(killed$biomass)) rep(NA_real_, n) else killed$biomass
    rowsum(cbind(cohorts = rep.int(1L, n), density = density), killed$species_code)
  })
}

.basefire_species_mortality <- function(log, killed_species, killed_density, area_mg) {
  rows <- lapply(which(log$cohorts_killed > 0L), function(i) {
    tally <- if (i <= length(killed_species)) killed_species[[i]]
    if (is.null(tally)) {
      return(data.frame(time = log$time[i], event_id = log$event_id[i],
                        species = NA_character_, cohorts_killed = log$cohorts_killed[i],
                        MgBiomassKilled = killed_density[i] * area_mg))
    }
    tally <- tally[order(rownames(tally)), , drop = FALSE]
    data.frame(time = log$time[i], event_id = log$event_id[i],
               species = rownames(tally), cohorts_killed = as.integer(tally[, "cohorts"]),
               MgBiomassKilled = unname(tally[, "density"]) * area_mg)
  })
  do.call(rbind, rows)
}

basefire_mortality <- function(result, extension = "fire", by_species = FALSE) {
  if (!isTRUE(by_species) && !identical(by_species, FALSE))
    .fire_stop("by_species must be TRUE or FALSE.")
  single_step <- inherits(result, "landisR_fire_basefire_landisR_result")
  if (!single_step && (!inherits(result, "landisR_coupled_run") ||
      !extension %in% names(result$events)))
    .fire_stop("Expected a basefire_step() result or a run with the named Base Fire extension.")
  simulation <- result$simulation
  # Biomass Succession stores cohort density in g/m2. Return Mg per fire.
  area_mg <- if (identical(simulation$growth_engine$id, "biomass52"))
    simulation$project$canonical$cell_length^2 / 1e6 else NA_real_

  if (single_step) {
    fire <- result$fire_result
    events <- list(list(log = fire$event_log,
                        killed_density = .basefire_killed_density(fire),
                        killed_species = if (by_species) .basefire_killed_species(fire)))
  } else {
    events <- result$events[[extension]]
  }
  logs <- lapply(events, function(event) {
    log <- event$log
    if (!is.data.frame(log) || !"cohorts_killed" %in% names(log))
      .fire_stop("The named extension did not return a Base Fire event log.")
    density <- event$killed_density
    if (is.null(density)) density <- ifelse(log$cohorts_killed == 0L, 0, NA_real_)
    if (by_species)
      return(.basefire_species_mortality(log, event$killed_species, density, area_mg))
    log$MgBiomassKilled <- density * area_mg
    log
  })
  if (by_species) {
    empty <- data.frame(time = integer(), event_id = integer(), species = character(),
                        cohorts_killed = integer(), MgBiomassKilled = numeric())
    logs <- c(list(empty), logs)
  } else if (!length(logs)) {
    return(data.frame(event_id = integer(), time = integer(), row = integer(),
      column = integer(), sites_checked = integer(), cohorts_killed = integer(),
      severity = numeric(), burned_sites = integer(), MgBiomassKilled = numeric()))
  }
  log <- do.call(rbind, unname(logs))
  rownames(log) <- NULL
  log
}
