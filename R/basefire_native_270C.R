# 270C : backend natif 265 permanent, code compile; charge apres les sources 270B.
# Aucun patch d'espace de noms a l'execution. Backend legacy par defaut
# tant que le nouveau paquet n'a pas prouve 18/18 + 24/24 contre 245.

.basefire_native_run_cpp_270C <- function(prepared_profile, current_time,
                                           cohort_env, cohort_sites, runtime,
                                           last_wind_year, rng, landis_rng = NULL) {
  timer <- .fire_timer("interface_cpp", current_time)
  # Hooks non standard : refus explicite plutot que changement silencieux du RNG.
  default_shuffle <- .basefire_default_shuffle
  if (!identical(rng$ignition_uniform,stats::runif) ||
      !identical(rng$generate_uniform,stats::runif) ||
      !identical(rng$next_double,stats::runif) ||
      !identical(rng$shuffle_fun,default_shuffle) ||
      !is.null(rng$wind_direction) || !is.null(rng$wind_speed)) {
    stop('270C: backend natif ne prend pas en charge ces hooks RNG/vent personnalises; utiliser backend legacy.',call.=FALSE)
  }
  pp <- prepared_profile$model$prepared
  spread <- pp$spread
  topo <- spread$site_topology
  mort <- pp$mortality
  sev <- pp$severity
  n <- length(spread$cell_id)
  if (n < 1L) stop('270C: paysage sans cellule')
  species_idx <- match(as.character(cohort_env$table$species_code),
                       as.character(mort$species_code))
  if (anyNA(species_idx)) stop('270C: espece de cohorte inconnue')
  sw <- if (is.null(last_wind_year)) numeric() else as.numeric(last_wind_year)
  has_wind <- isTRUE(sev$has_wind_state) || !is.null(last_wind_year)
  size_one <- function(region_index) {
    .basefire_call_size(
      region_index=region_index,
      size_fun=prepared_profile$model$size,
      size_prepared=pp$size,
      exponential_draw=rng$exponential_draw
    )
  }
  v7_config <- if (!is.null(landis_rng)) {
    scan_order <- topo$scan_order_v7
    if (is.null(scan_order))
      scan_order <- order(prepared_profile$topology$row, prepared_profile$topology$col, method = "radix")
    c(landis_rng, list(scan_order = scan_order, mean_size = pp$size$mean_size,
      min_size = pp$size$min_size, max_size = pp$size$max_size, cell_area = pp$size$cell_area,
      timing = .fire_timing_enabled()))
  }
  timer("arguments")
  result <- basefire_timestep_kernel_260(
    as.integer(current_time), as.numeric(prepared_profile$profile$timestep),
    as.numeric(pp$ignition$ignition_probability),
    as.numeric(pp$ignition$fire_spread_age),
    as.integer(spread$cell_id), as.integer(spread$region_index),
    as.integer(spread$n_regions), topo$neighbors_full,
    topo$second_neighbors_full, as.integer(topo$full_to_active),
    as.logical(topo$full_active), as.numeric(pp$burn_probability$fire_spread_age),
    sev$fuel_curve, sev$wind_curve, has_wind, sw,
    as.numeric(runtime$last_fire_year), as.logical(runtime$event_occupied),
    as.numeric(runtime$last_severity), cohort_sites, cohort_env,
    as.numeric(cohort_env$table$age), as.integer(species_idx),
    as.integer(mort$species_fire_tolerance),
    as.numeric(mort$species_longevity),
    as.numeric(mort$damage_max_age),
    as.integer(mort$damage_severity_tolerance_diff),
    size_one, v7_config
  )
  timer("noyau")
  if (!is.null(result$timing))
    .fire_add_timing("cpp", current_time, names(result$timing), result$timing)
  result
}

.basefire_run_timestep_native_270C <- function(
    prepared_profile,
    current_time,
    cohorts,
    rng = list(),
    landis_rng = NULL) {

  timer <- .fire_timer("moteur", current_time)
  .basefire_validate_prepared_profile(
    prepared_profile
  )

  if (!is.numeric(current_time) ||
      length(current_time) != 1L ||
      is.na(current_time) ||
      !is.finite(current_time)) {
    .fire_stop(
      'current_time must be one finite numeric value.'
    )
  }

  current_time <- as.integer(current_time)

  timer("validation")
  prepared_profile <-
    .basefire_apply_dynamic_region_update(
      prepared_profile,
      current_time
    )

  timer("carte_dynamique")
  topology <- prepared_profile$topology

  .basefire_validate_cohort_table(
    cohorts,
    topology
  )

  rng <- .basefire_normalize_rng(rng)

  timer("validation_cohortes")
  cohort_env <- new.env(parent = emptyenv())
  cohort_env$table <- cohorts
  rownames(cohort_env$table) <- NULL
  # Masque de mortalite stable pendant toute la decennie.
  cohort_env$alive <- rep.int(TRUE, nrow(cohort_env$table))
  cohort_ids <- cohorts$cell_id
  # Preserve the native path's original grouping for nonstandard cell IDs.
  cohort_sites <- if (is.integer(cohort_ids) && !is.object(cohort_ids)) {
    match(cohort_ids, topology$cell_id, nomatch = 0L)
  } else {
    match(as.character(cohort_ids), as.character(topology$cell_id), nomatch = 0L)
  }

  timer("index_cohortes")
  n <- length(topology$cell_id)

  last_fire_year <- prepared_profile$model$prepared$
    ignition$last_fire_year

  last_severity <- rep(NA_real_, n)

  if (!is.null(prepared_profile$fire_state)) {
    state_idx <- match(
      topology$cell_id,
      prepared_profile$fire_state$cells$cell_id
    )

    if (anyNA(state_idx)) {
      .fire_stop(
        'Prepared Base Fire state is missing active cells.'
      )
    }

    last_severity <-
      prepared_profile$fire_state$cells$
        last_severity[state_idx]
  }

  runtime <- new.env(parent = emptyenv())
  runtime$last_fire_year <- as.numeric(last_fire_year)
  runtime$last_severity <- as.numeric(last_severity)
  runtime$event_occupied <- rep(FALSE, n)
  runtime$event_owner <- integer(n)
  runtime$severity <- numeric(n)
  runtime$disturbed <- rep(FALSE, n)

  last_wind_year <- prepared_profile$model$prepared$
    severity$last_wind_year

  severity_fun <- function(cell_index, current_time) {
    prepared_profile$model$severity(
      cell_index = cell_index,
      current_time = current_time,
      prepared = prepared_profile$model$prepared$severity,
      last_fire_year = runtime$last_fire_year,
      last_wind_year = last_wind_year,
      diagnostics = FALSE
    )
  }

  damage_fun <- function(
      cell_index,
      severity,
      current_time) {

    cell_id <- topology$cell_id[cell_index]

    cohort_rows <- which(
      cohort_env$table$cell_id == cell_id
    )

    if (!length(cohort_rows)) {
      return(0L)
    }

    local <- cohort_env$table[
      cohort_rows,
      ,
      drop = FALSE
    ]

    mortality <- .basefire_apply_mortality_component(
      cohorts = local,
      severity = severity,
      mortality_fun = prepared_profile$model$mortality,
      mortality_prepared =
        prepared_profile$model$prepared$mortality
    )

    if (any(mortality$kill)) {
      global_kill <- cohort_rows[mortality$kill]
      cohort_env$table <- cohort_env$table[
        -global_kill,
        ,
        drop = FALSE
      ]
      rownames(cohort_env$table) <- NULL
    }

    mortality$cohorts_killed
  }

  raw_events <- list()
  events <- list()
  log_rows <- list()

  summary_region_sites <- integer(
    prepared_profile$parameter_cache$n_regions
  )

  summary_total_sites <- 0L
  summary_event_count <- 0L

  # 270C: delegation explicite au noyau 265 compile dans ce paquet.
  timer("preparation_etat")
  fused270C <- .basefire_native_run_cpp_270C(
    prepared_profile=prepared_profile, current_time=current_time,
    cohort_env=cohort_env, cohort_sites=cohort_sites, runtime=runtime,
    last_wind_year=last_wind_year, rng=rng, landis_rng=landis_rng
  )
  timer("interface_et_noyau")
  runtime$last_fire_year <- fused270C$last_fire_year
  runtime$event_occupied <- fused270C$event_occupied
  runtime$last_severity <- fused270C$last_severity
  runtime$event_owner <- fused270C$event_owner
  runtime$severity <- fused270C$severity
  runtime$disturbed <- fused270C$disturbed
  raw_events <- fused270C$raw_events
  summary_event_count <- length(raw_events)
  for (event_id in seq_along(raw_events)) {
    spread_result <- raw_events[[event_id]]
    if (any(c('event_member','event_occupied','last_fire_year') %in%
            names(spread_result))) {
      stop('270C: instantane complet detecte dans raw_events')
    }
    kill_idx <- spread_result$killed_cohort_rows
    if (length(kill_idx) != spread_result$cohorts_killed ||
        anyDuplicated(kill_idx) || anyNA(kill_idx) ||
        any(kill_idx < 1L | kill_idx > nrow(cohort_env$table))) {
      stop('270C: indices des cohortes mortes incoherents')
    }
    affected <- cohort_env$table[kill_idx,,drop=FALSE]
    rownames(affected) <- NULL
    affected$cohort_row_before_fire <- as.integer(kill_idx)
    affected$effect <- rep.int('killed',length(kill_idx))
    spread_result$affected_cohorts <- affected
    raw_events[[event_id]] <- spread_result
    standard_event <- .basefire_spread_to_fire_event(
      spread_result=spread_result,prepared_profile=prepared_profile,
      current_time=current_time,event_id=event_id
    )
    if (!is.null(standard_event)) events[[length(events)+1L]] <- standard_event
    log_row <- .basefire_event_log_row(
      spread_result=spread_result, prepared_profile=prepared_profile,
      current_time=current_time, event_id=event_id
    )
    if (!is.null(log_row)) {
      log_rows[[length(log_rows)+1L]] <- log_row
      summary_total_sites <- summary_total_sites + spread_result$total_sites_in_event
      summary_region_sites <- summary_region_sites + spread_result$sites_in_event
    }
  }

  if (length(log_rows)) {
    event_log <- do.call(rbind, log_rows)
    rownames(event_log) <- NULL
  } else {
    event_log <- data.frame(
      event_id = integer(),
      time = integer(),
      row = integer(),
      column = integer(),
      sites_checked = integer(),
      cohorts_killed = integer(),
      severity = numeric(),
      burned_sites = integer(),
      stringsAsFactors = FALSE
    )
  }

  names(summary_region_sites) <-
    prepared_profile$parameter_cache$region_name

  summary <- list(
    time = current_time,
    total_sites_burned = as.integer(
      summary_total_sites
    ),
    number_fires = as.integer(
      summary_event_count
    ),
    region_sites = summary_region_sites
  )

  timer("journaux")
  state <- fire_state(
    cell_id = topology$cell_id,
    last_fire_year = runtime$last_fire_year,
    last_severity = runtime$last_severity
  )

  if (!is.null(last_wind_year)) {
    state$cells$last_wind_year <- last_wind_year
  }

  # Carry the updated state into the prepared profile so the same object can
  # be used for the next timestep without reinitializing fire history.
  prepared_profile$fire_state <- state

  prepared_profile$model$prepared$ignition$
    last_fire_year <- runtime$last_fire_year

  prepared_profile$model$prepared$burn_probability$
    last_fire_year <- runtime$last_fire_year

  prepared_profile$model$prepared$severity$
    last_fire_year <- runtime$last_fire_year

  # Compaction differee apres TOUS les feux du timestep.
  cohort_env$table <- cohort_env$table[cohort_env$alive,,drop=FALSE]
  rownames(cohort_env$table) <- NULL

  out <- list(
    time = current_time,
    events = events,
    raw_events = raw_events,
    event_log = event_log,
    summary = summary,
    cohorts = cohort_env$table,
    fire_state = state,
    severity = data.frame(
      cell_id = topology$cell_id,
      severity = runtime$severity,
      disturbed = runtime$disturbed,
      event_owner = runtime$event_owner
    ),
    prepared_profile = prepared_profile
  )

  if (!is.null(landis_rng)) out$landis_rng <- fused270C$landis_rng
  timer("etat_et_compaction")
  class(out) <- 'landisR_fire_basefire_timestep_result'
  out
}


# Signature historique stricte : option seulement, sans nouvel argument.
.basefire_run_timestep <- function(prepared_profile, current_time,
                                   cohorts, rng = list(), landis_rng = NULL) {
  backend <- getOption("landisR.fire.backend", NULL)
  if (is.null(backend)) {
    rr <- .basefire_normalize_rng(rng)
    compatible <- identical(prepared_profile$model$name, "basefire-strict") &&
      identical(rr$ignition_uniform, stats::runif) &&
      identical(rr$generate_uniform, stats::runif) &&
      identical(rr$next_double, stats::runif) &&
      identical(rr$shuffle_fun, .basefire_default_shuffle) &&
      is.null(rr$wind_direction) && is.null(rr$wind_speed)
    backend <- if (compatible) "native" else "legacy"
  }
  if (!identical(backend,'legacy') && !identical(backend,'native'))
    stop('landisR.fire.backend doit etre legacy ou native.',call.=FALSE)
  if (!is.null(landis_rng) && backend != "native")
    stop("The LANDIS-II v7 fire stream requires the native backend.", call. = FALSE)
  if (identical(backend,'native')) {
    return(.basefire_run_timestep_native_270C(prepared_profile,
               current_time, cohorts, rng, landis_rng))
  }
  .basefire_run_timestep_legacy_270C(prepared_profile, current_time,
                                     cohorts, rng)
}

.basefire_native_status_270B <- function() {
  ns <- asNamespace('landisR.fire')
  expected <- c('basefire_valid_cohorts_fast_263',
                'basefire_timestep_kernel_260',
                'basefire_event_kernel_247',
                'basefire_rng_sequence_249',
                'basefire_kernel_stats_249')
  available <- vapply(expected,exists,logical(1),envir=ns,
                      mode='function', inherits=FALSE)
  backend <- getOption("landisR.fire.backend", "native")
  list(package_version=as.character(utils::packageVersion('landisR.fire')),
       compiled_entrypoints=available,
       native_validation_default=isTRUE(getOption('landisR.fire.native_validation',TRUE)),
       fused_native_backend_active=identical(backend,'native'),
       production_fire_backend=backend,
       information='Native default for strict Base Fire with standard RNG; legacy available for custom hooks.')
}
