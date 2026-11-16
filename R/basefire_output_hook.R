# Public Base Fire convenience API with strictly optional writing.
# This is a compiled-package source extension, not a runtime namespace patch.
.basefire_step_without_outputs_271 <- basefire_step

#' @inheritParams basefire_write_outputs
#' @export
basefire_step <- function(simulation, prepared, rng = list(),
                          stochastic_serotiny = TRUE, metadata = list(),
                          progress = FALSE, verbose = FALSE,
                          output_dir = NULL, output_formats = c('csv','rds'),
                          raster_template = NULL, cell_length = NULL,
                          overwrite_outputs = FALSE, year = simulation$year) {
  out <- .basefire_step_without_outputs_271(
    simulation=simulation, prepared=prepared, rng=rng,
    stochastic_serotiny=stochastic_serotiny, metadata=metadata,
    progress=progress, verbose=verbose, year=year
  )
  if (!is.null(output_dir)) {
    out$output_files <- basefire_write_outputs(
      result=out, output_dir=output_dir, formats=output_formats,
      raster_template=raster_template, cell_length=cell_length,
      overwrite=overwrite_outputs
    )
    if (isTRUE(verbose)) cat('[landisR.fire] Outputs written to ',output_dir,'\n',sep='')
  }
  out
}
