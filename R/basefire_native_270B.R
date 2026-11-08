# 270B: Code permanent du paquet de test. Pas de monkey patch dans une session.
# Le backend de propagation natif 265 est compile mais inactif par defaut.
# Le validateur 263B (strictement conservateur) peut accelerer le controle.
.basefire_validate_cohort_table <- function(cohorts, topology) {
  if (isTRUE(getOption('landisR.fire.native_validation', TRUE))) {
    fast <- tryCatch(basefire_valid_cohorts_fast_263(cohorts, topology$cell_id),
                     error = function(e) FALSE)
    if (isTRUE(fast)) return(invisible(TRUE))
  }
  .basefire_validate_cohort_table_legacy_270B(cohorts, topology)
}

.basefire_native_status_270B <- function() {
  expected <- c('basefire_valid_cohorts_fast_263',
                'basefire_timestep_kernel_260',
                'basefire_event_kernel_247',
                'basefire_rng_sequence_249',
                'basefire_kernel_stats_249')
  ns <- asNamespace('landisR.fire')
  available <- vapply(expected, exists, logical(1), envir = ns,
                      mode = 'function', inherits = FALSE)
  list(package_version = as.character(utils::packageVersion('landisR.fire')),
       compiled_entrypoints = available,
       native_validation_default = isTRUE(getOption('landisR.fire.native_validation', TRUE)),
       fused_native_backend_active = FALSE,
       production_fire_backend = 'legacy R (safe default)',
       information = 'Rcpp compiled into package; fused 265 requires separate parity tests before activation')
}

# API compatible: nouveaux arguments ajoutes apres les anciens.
# La progression marque les phases de l'appel public, pas les cellules individuelles.
basefire_step <- function(simulation, prepared, rng = list(),
                          stochastic_serotiny = TRUE, metadata = list(),
                          progress = FALSE, verbose = FALSE) {
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress))
    stop('progress must be TRUE or FALSE.', call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop('verbose must be TRUE or FALSE.', call. = FALSE)
  if (!progress && !verbose) {
    return(.basefire_step_landisR(simulation = simulation,
       prepared_profile = prepared, rng = rng,
       stochastic_serotiny = stochastic_serotiny, metadata = metadata))
  }
  # Consommation RNG nulle: aucune RNG ne doit intervenir dans l'interface.
  bar <- if (progress) utils::txtProgressBar(min=0, max=2, style=3) else NULL
  if (!is.null(bar)) on.exit(close(bar), add=TRUE)
  if (verbose) cat('[landisR.fire] Debut du pas Base Fire\n')
  if (!is.null(bar)) utils::setTxtProgressBar(bar, 1)
  ans <- .basefire_step_landisR(simulation = simulation,
    prepared_profile = prepared, rng = rng,
    stochastic_serotiny = stochastic_serotiny, metadata = metadata)
  if (!is.null(bar)) utils::setTxtProgressBar(bar, 2)
  if (verbose) cat('[landisR.fire] Fin du pas Base Fire\n')
  ans
}
