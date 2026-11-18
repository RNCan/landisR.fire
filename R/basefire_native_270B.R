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
