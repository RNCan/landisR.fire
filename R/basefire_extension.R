# Keep the R draws for the fire engine in its own checkpointed stream.
.basefire_coupled_rng <- function(rng_state, step) {
  previous <- .GlobalEnv$.Random.seed
  on.exit({
    if (is.null(previous)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", previous, envir = .GlobalEnv)
    }
  })
  assign(".Random.seed", rng_state, envir = .GlobalEnv)
  result <- step()
  list(result = result, rng_state = .GlobalEnv$.Random.seed)
}

basefire_extension <- function(parameters, fire_region_map_code = NULL, cell_area = NULL,
                               fire_state = NULL, dynamic_region_map_reader = NULL,
                               time_of_last_wind = NULL, parameter_schedule = NULL,
                               seed = NULL, stochastic_serotiny = TRUE,
                               rng_mode = c("independent", "landis_v7")) {
  rng_mode <- match.arg(rng_mode)
  shared <- rng_mode == "landis_v7"
  if (shared && !is.null(seed))
    .fire_stop("Use the simulation seed with the shared LANDIS-II v7 fire stream.")
  force(parameters)
  force(fire_region_map_code)
  force(cell_area)
  force(fire_state)
  force(dynamic_region_map_reader)
  force(time_of_last_wind)
  force(parameter_schedule)
  force(seed)
  force(stochastic_serotiny)

  initialize <- function(simulation) {
    if (shared && is.null(simulation$landis_rng))
      .fire_stop("The shared fire stream requires recruitment_rng = 'landis_v7'.")
    area <- if (is.null(cell_area)) simulation$project$canonical$cell_length^2 / 10000 else cell_area
    prepared <- .basefire_prepare_landisR(
      simulation, parameters, fire_region_map_code, area,
      fire_state, dynamic_region_map_reader, time_of_last_wind,
      shared_geometry = shared
    )
    rng_state <- if (!shared) simulation$rng_state
    if (!is.null(seed)) {
      if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
          seed < 0 || seed > .Machine$integer.max || seed != trunc(seed))
        .fire_stop("Base Fire seed must be a nonnegative integer.")
      rng_state <- .basefire_coupled_rng(rng_state, function() set.seed(seed))$rng_state
    }
    list(prepared = prepared, rng_state = rng_state)
  }

  step <- function(simulation, state, year) {
    report <- getOption("landisR.progress.report")
    prepared <- state$prepared
    if (!is.null(parameter_schedule)) {
      if (!is.null(report)) report(0, "Updating fire parameters")
      prepared <- basefire_update_parameters(prepared, parameter_schedule, current_time = year)
    }
    if (shared) {
      result <- .basefire_step_landisR(simulation, prepared,
        stochastic_serotiny = stochastic_serotiny, year = year, use_landis_rng = TRUE)
    } else {
      draw <- .basefire_coupled_rng(state$rng_state, function() {
        basefire_step(simulation, prepared, stochastic_serotiny = stochastic_serotiny,
                      progress = FALSE, year = year)
      })
      result <- draw$result
      state$rng_state <- draw$rng_state
    }
    if (!is.null(report)) report(0.97, "Collecting fire results")
    state$prepared <- result$prepared_profile
    severity <- result$fire_result$severity
    burned <- severity[severity$disturbed & severity$severity > 0,
                       c("cell_id", "severity"), drop = FALSE]
    rownames(burned) <- NULL
    list(simulation = result$simulation,
         state = state,
         events = list(burned = burned, log = result$fire_result$event_log,
                       summary = result$fire_result$summary,
                       killed_density = .basefire_killed_density(result$fire_result),
                       killed_species = .basefire_killed_species(result$fire_result)))
  }

  .basefire_require_landisR()
  structure(landisR::landis_extension(step, initialize, type = "disturbance"),
            progress_label = "Base Fire")
}
