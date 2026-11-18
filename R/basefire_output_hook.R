# Public Base Fire convenience API with strictly optional writing.
# This is a compiled-package source extension, not a runtime namespace patch.
#' @inheritParams basefire_write_outputs
#' @export
basefire_step <- function(simulation, prepared, rng = list(),
                          stochastic_serotiny = TRUE, metadata = list(),
                          progress = getOption("landisR.progress", interactive()), verbose = FALSE,
                          output_dir = NULL, output_formats = c('csv','rds'),
                          raster_template = NULL, cell_length = NULL,
                          overwrite_outputs = FALSE, year = simulation$year) {
  if (!identical(progress, TRUE) && !identical(progress, FALSE))
    stop("progress must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  if (verbose) cat("[landisR.fire] Debut du pas Base Fire\n")
  report <- getOption("landisR.progress.report")
  bar <- landisR::landis_progress_bar(progress && is.null(report), "Base Fire")
  if (!is.null(bar)) {
    on.exit(bar$close(), add = TRUE)
    report <- function(fraction, label)
      bar$update(fraction, paste("Year", year, "|", label))
  }
  out <- .basefire_step_landisR(
    simulation = simulation, prepared_profile = prepared, rng = rng,
    stochastic_serotiny = stochastic_serotiny, metadata = metadata,
    year = year, report_progress = report
  )
  if (!is.null(output_dir)) {
    if (!is.null(report)) report(0.98, "Writing outputs")
    out$output_files <- basefire_write_outputs(
      result=out, output_dir=output_dir, formats=output_formats,
      raster_template=raster_template, cell_length=cell_length,
      overwrite=overwrite_outputs
    )
  }
  if (!is.null(bar)) {
    report(1, "Complete")
    bar$close()
  }
  if (verbose) {
    cat("[landisR.fire] Fin du pas Base Fire\n")
    if (!is.null(output_dir))
      cat("[landisR.fire] Outputs written to ", output_dir, "\n", sep = "")
  }
  out
}
