.basefire_same_raster_grid <- function(reference, raster) {
  if (terra::compareGeom(reference, raster, res = TRUE, stopOnError = FALSE))
    return(TRUE)
  if (!terra::compareGeom(reference, raster, crs = FALSE, res = TRUE,
                          stopOnError = FALSE))
    return(FALSE)

  # The 6est Lambert Conformal Conic maps list the same standard parallels
  # in opposite order. Their order does not affect the projection.
  proj_tokens <- function(x) {
    strsplit(trimws(terra::crs(x, proj = TRUE)), '[[:space:]]+')[[1L]]
  }
  ref_proj <- proj_tokens(reference)
  map_proj <- proj_tokens(raster)
  if (!('+proj=lcc' %in% ref_proj && '+proj=lcc' %in% map_proj))
    return(FALSE)
  lcc_key <- function(x) sort(sub('^[+]lat_[12]=', '+lat_parallel=', x))
  identical(lcc_key(ref_proj), lcc_key(map_proj))
}

.basefire_raster_map_reader <- function(spec, nrow, ncol) {
  force(spec)
  force(nrow)
  force(ncol)

  function(map_name, current_time = NULL, landscape) {
    if (!requireNamespace('terra', quietly = TRUE))
      .fire_stop('Install terra to read Base Fire region rasters.')
    path <- .basefire_resolve_path(chartr('\\', '/', map_name), spec$directory)
    if (!file.exists(path))
      .fire_stop('Base Fire region raster not found: ', path)
    raster <- terra::rast(path)
    if (terra::nlyr(raster) != 1L || terra::nrow(raster) != nrow ||
        terra::ncol(raster) != ncol)
      .fire_stop('Expected a one-layer raster matching the landscape dimensions: ', path)

    if (!identical(path, spec$template)) {
      if (!file.exists(spec$template))
        .fire_stop('Base Fire reference raster not found: ', spec$template)
      if (!.basefire_same_raster_grid(terra::rast(spec$template), raster))
        .fire_stop('Base Fire region raster is not aligned to the reference grid: ', path)
    }

    # terra enumerates cells row by row. Use canonical rows and columns.
    index <- (landscape$row - 1L) * ncol + landscape$col
    codes <- terra::values(raster, mat = FALSE)[index]
    codes[is.na(codes)] <- 0
    .basefire_integer_like(codes, 'Base Fire raster region codes')
  }
}
