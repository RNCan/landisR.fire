# Base Fire output writer, 271. All writing is deterministic and opt-in.
# Does not modify forest state, fire events, random numbers or Rcpp kernels.

.basefire_out_stop <- function(...) stop(..., call. = FALSE)

.basefire_out_result <- function(x) {
  if (!is.list(x)) .basefire_out_stop('result must be a Base Fire result list.')
  fr <- if (!is.null(x$fire_result)) x$fire_result else x
  if (!is.list(fr) || !is.data.frame(fr$event_log) ||
      !is.data.frame(fr$severity) || !is.list(fr$summary) ||
      !is.list(fr$fire_state) || !is.list(fr$events)) {
    .basefire_out_stop('Unsupported Base Fire result: expected event_log, severity, summary, fire_state and events.')
  }
  prep <- if (!is.null(x$prepared_profile)) x$prepared_profile else fr$prepared_profile
  list(fire = fr, prepared = prep)
}

.basefire_out_year <- function(fr) {
  y <- if (!is.null(fr$time)) fr$time else fr$summary$time
  if (!is.numeric(y) || length(y) != 1L || is.na(y) ||
      !is.finite(y) || abs(y - round(y)) > 1e-8) {
    .basefire_out_stop('The Base Fire result must have one integer simulation year.')
  }
  as.integer(y)
}

.basefire_out_event_cells <- function(events, year) {
  ans <- vector('list', length(events))
  k <- 0L
  for (i in seq_along(events)) {
    ev <- events[[i]]
    if (!is.list(ev) || is.null(ev$severity) || !is.data.frame(ev$severity) ||
        !all(c('cell_id','severity') %in% names(ev$severity))) next
    ec <- ev$severity[, c('cell_id','severity'), drop = FALSE]
    if (!nrow(ec)) next
    id <- if (!is.null(ev$metadata$event_id)) ev$metadata$event_id else i
    k <- k + 1L
    ans[[k]] <- data.frame(year = rep.int(year, nrow(ec)),
                           event_id = rep.int(as.integer(id), nrow(ec)),
                           cell_id = as.integer(ec$cell_id),
                           severity = as.numeric(ec$severity))
  }
  if (!k) return(data.frame(year = integer(), event_id = integer(),
                            cell_id = integer(), severity = numeric()))
  out <- do.call(rbind, ans[seq_len(k)])
  rownames(out) <- NULL
  out
}

.basefire_out_region <- function(fr, cell_area) {
  regional <- fr$summary$region_sites
  if (is.null(regional)) regional <- numeric()
  values <- as.numeric(regional)
  nms <- names(regional)
  if (is.null(nms)) nms <- as.character(seq_along(values))
  data.frame(fire_region = nms,
             burned_sites = values,
             burned_area_cell_area_units = if (is.finite(cell_area)) values * cell_area else rep(NA_real_, length(values)),
             stringsAsFactors = FALSE)
}

.basefire_out_scalar <- function(fr, cell_area) {
  s <- fr$summary
  n_fires <- if (!is.null(s$number_fires)) as.integer(s$number_fires) else length(fr$raw_events)
  n_sites <- if (!is.null(s$total_sites_burned)) as.numeric(s$total_sites_burned) else NA_real_
  data.frame(year = .basefire_out_year(fr),
             number_fires = n_fires,
             total_sites_burned = n_sites,
             total_burned_area_cell_area_units = if (is.finite(cell_area)) n_sites * cell_area else NA_real_)
}

.basefire_out_grid <- function(fr, prepared) {
  if (is.null(prepared) || !is.data.frame(prepared$landscape) ||
      !all(c('cell_id','row','col') %in% names(prepared$landscape))) {
    .basefire_out_stop('Raster output requires prepared_profile$landscape with cell_id, row and col.')
  }
  land <- prepared$landscape
  ids <- as.integer(land$cell_id)
  rr <- as.integer(land$row)
  cc <- as.integer(land$col)
  if (!length(ids) || anyNA(c(ids,rr,cc)) || anyDuplicated(ids) ||
      anyDuplicated(data.frame(row=rr,col=cc))) {
    .basefire_out_stop('Invalid landscape cell identifiers or grid coordinates.')
  }
  if (!is.data.frame(fr$severity) || !all(c('cell_id','severity','disturbed') %in% names(fr$severity))) {
    .basefire_out_stop('Fire result lacks the cell severity table.')
  }
  sev <- fr$severity
  if (anyDuplicated(sev$cell_id)) .basefire_out_stop('Duplicate cell IDs in severity table.')
  ix <- match(as.integer(sev$cell_id), ids)
  if (anyNA(ix)) .basefire_out_stop('Fire severity contains cells absent from the landscape.')
  active <- if ('active' %in% names(land)) as.logical(land$active) else rep(TRUE,length(ids))
  if (anyNA(active)) .basefire_out_stop('Invalid active-cell mask.')
  if (any(!active[ix])) .basefire_out_stop('Severity contains inactive cells.')
  y <- .basefire_out_year(fr)
  hist <- fr$fire_state$cells
  if (!is.data.frame(hist) || !all(c('cell_id','last_fire_year') %in% names(hist))) {
    .basefire_out_stop('Fire result lacks fire_state$cells history data.')
  }
  hi <- match(ids, hist$cell_id)
  years <- as.numeric(hist$last_fire_year[hi])
  years[!active] <- NA_real_
  elapsed <- y - years
  elapsed[is.na(years)] <- NA_real_
  if (any(elapsed[!is.na(elapsed)] < 0)) .basefire_out_stop('Future last_fire_year encountered.')
  nrows <- max(rr) - min(rr) + 1L
  ncols <- max(cc) - min(cc) + 1L
  rows <- rr - min(rr) + 1L
  cols <- cc - min(cc) + 1L
  make_grid <- function(x, background = NA_real_) {
    m <- matrix(background, nrow=nrows, ncol=ncols)
    m[cbind(rows,cols)] <- x
    m
  }
  sev_values <- rep(NA_real_, length(ids))
  sev_values[active] <- 0
  sev_values[ix] <- as.numeric(sev$severity)
  ids_values <- rep(NA_real_, length(ids))
  ids_values[active] <- 0
  if ('event_owner' %in% names(sev)) {
    ids_values[ix] <- as.numeric(sev$event_owner)
  } else if (length(fr$events)) {
    .basefire_out_stop('Event-ID raster needs severity$event_owner.')
  }
  list(fire_severity = make_grid(sev_values),
       fire_event_id = make_grid(ids_values),
       time_since_fire = make_grid(elapsed))
}

.basefire_out_template <- function(template, cell_length, nrows, ncols) {
  if (!requireNamespace('terra', quietly = TRUE))
    .basefire_out_stop('GeoTIFF export requires the optional terra package.')
  if (!is.null(template)) {
    t <- if (is.character(template) && length(template) == 1L) terra::rast(template) else template
    if (!inherits(t, 'SpatRaster')) .basefire_out_stop('raster_template must be a SpatRaster or a raster path.')
    if (terra::nlyr(t) != 1L || terra::nrow(t) != nrows || terra::ncol(t) != ncols)
      .basefire_out_stop('raster_template dimensions must match the full landscape grid.')
    return(list(template = t, relative = FALSE))
  }
  if (!is.numeric(cell_length) || length(cell_length) != 1L ||
      !is.finite(cell_length) || cell_length <= 0)
    .basefire_out_stop('GeoTIFF requires raster_template for georeferencing, or cell_length for relative grid coordinates.')
  r <- terra::rast(nrows=nrows,ncols=ncols,xmin=0,xmax=ncols*cell_length,
                   ymin=0,ymax=nrows*cell_length)
  list(template = r, relative = TRUE)
}

#' Write Base Fire outputs without changing the simulation or RNG
#'
#' @param result Return value of basefire_step(), or its fire_result member.
#' @param output_dir Parent folder for annual outputs.
#' @param formats Output formats, a subset of c('csv','rds','tif').
#' @param raster_template Optional one-layer terra SpatRaster or file path with
#'   correct grid dimensions, extent and CRS. Required for true georeferencing.
#' @param cell_length Optional cell size for *relative*, non-georeferenced rasters.
#' @param overwrite Allow replacement of an existing year's output folder.
#' @return Invisibly, a manifest of files written.
#' @export
basefire_write_outputs <- function(result, output_dir,
                                   formats = c('csv','rds'),
                                   raster_template = NULL,
                                   cell_length = NULL,
                                   overwrite = FALSE) {
  parsed <- .basefire_out_result(result)
  fr <- parsed$fire
  prep <- parsed$prepared
  year <- .basefire_out_year(fr)
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir))
    .basefire_out_stop('output_dir must be a nonempty path.')
  if (!is.character(formats) || !length(formats) || anyNA(formats) ||
      anyDuplicated(formats) || !all(formats %in% c('csv','rds','tif')))
    .basefire_out_stop("formats must contain distinct values among 'csv', 'rds', 'tif'.")
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite))
    .basefire_out_stop('overwrite must be TRUE or FALSE.')
  grids <- NULL
  r_spec <- NULL
  if ('tif' %in% formats) {
    grids <- .basefire_out_grid(fr, prep)
    r_spec <- .basefire_out_template(raster_template, cell_length,
                                   nrow(grids[[1L]]), ncol(grids[[1L]]))
  }
  cell_area <- NA_real_
  if (!is.null(prep$model$prepared$size$cell_area))
    cell_area <- as.numeric(prep$model$prepared$size$cell_area)[1L]
  if (!length(cell_area) || !is.finite(cell_area)) cell_area <- NA_real_
  output_dir <- path.expand(output_dir)
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  if (!dir.exists(output_dir)) .basefire_out_stop('Cannot create output_dir: ', output_dir)
  year_dir <- file.path(output_dir, sprintf('year_%04d', year))
  if (file.exists(year_dir) && !overwrite)
    .basefire_out_stop('Output folder already exists; use overwrite=TRUE: ', year_dir)
  tmp <- file.path(output_dir, paste0('.year_',sprintf('%04d',year),'.tmp_',Sys.getpid()))
  if (file.exists(tmp)) .basefire_out_stop('Temporary output folder exists: ',tmp)
  if (!dir.create(tmp)) .basefire_out_stop('Unable to create temporary output folder.')
  backed_up <- NULL
  published <- FALSE
  on.exit({
    if (dir.exists(tmp)) unlink(tmp, recursive=TRUE, force=TRUE)
    if (!published && !is.null(backed_up) && dir.exists(backed_up) &&
        !dir.exists(year_dir)) file.rename(backed_up,year_dir)
  }, add=TRUE)
  files <- character()
  categories <- character()
  register <- function(name, category) {
    files <<- c(files,name)
    categories <<- c(categories,category)
  }
  if ('csv' %in% formats) {
    ev <- fr$event_log
    ev$year <- rep.int(year,nrow(ev))
    utils::write.csv(ev,file.path(tmp,'fire_events.csv'),row.names=FALSE,na='')
    register('fire_events.csv','events')
    utils::write.csv(.basefire_out_event_cells(fr$events,year),
                     file.path(tmp,'fire_event_cells.csv'),row.names=FALSE,na='')
    register('fire_event_cells.csv','event_cells')
    utils::write.csv(.basefire_out_scalar(fr,cell_area),
                     file.path(tmp,'fire_summary.csv'),row.names=FALSE,na='')
    register('fire_summary.csv','summary')
    utils::write.csv(.basefire_out_region(fr,cell_area),
                     file.path(tmp,'fire_regions.csv'),row.names=FALSE,na='')
    register('fire_regions.csv','regions')
    utils::write.csv(fr$fire_state$cells,file.path(tmp,'fire_history.csv'),
                     row.names=FALSE,na='')
    register('fire_history.csv','history')
  }
  if ('rds' %in% formats) {
    saveRDS(fr$raw_events,file.path(tmp,'fire_raw_events.rds'),version=3)
    register('fire_raw_events.rds','sparse_events')
    saveRDS(fr$fire_state,file.path(tmp,'fire_state.rds'),version=3)
    register('fire_state.rds','fire_state')
  }
  if ('tif' %in% formats) {
    dir.create(file.path(tmp,'rasters'))
    for (nm in names(grids)) {
      r <- terra::setValues(terra::rast(r_spec$template), as.vector(t(grids[[nm]])))
      names(r) <- nm
      f <- file.path('rasters',paste0(nm,'.tif'))
      terra::writeRaster(r,file.path(tmp,f),overwrite=FALSE)
      register(f,'raster')
    }
    if (r_spec$relative) {
      writeLines(c('These rasters have RELATIVE grid coordinates, not a verified geographic CRS.',
                   'Supply raster_template for georeferenced GeoTIFF exports.',
                   'Grid row increases downward; row 1 is top.'),
                 file.path(tmp,'RASTERS_RELATIVE_COORDINATES.txt'))
      register('RASTERS_RELATIVE_COORDINATES.txt','metadata')
    }
  }
  manifest <- data.frame(category=categories,file=files,
                         size_bytes=as.numeric(file.info(file.path(tmp,files))$size),
                         stringsAsFactors=FALSE)
  utils::write.csv(manifest,file.path(tmp,'manifest.csv'),row.names=FALSE)
  if (dir.exists(year_dir)) {
    backed_up <- file.path(output_dir, paste0('.year_',sprintf('%04d',year),
                                              '.backup_',Sys.getpid()))
    if (file.exists(backed_up)) .basefire_out_stop('Existing backup folder: ',backed_up)
    if (!file.rename(year_dir,backed_up)) .basefire_out_stop('Cannot back up existing output directory.')
  }
  if (!file.rename(tmp,year_dir)) .basefire_out_stop('Failed to publish complete output directory.')
  published <- TRUE
  if (!is.null(backed_up) && dir.exists(backed_up)) unlink(backed_up,recursive=TRUE)
  manifest$path <- file.path(year_dir,manifest$file)
  invisible(manifest)
}
