# Test fixture for package checks and 271 smoke; no fire simulation or RNG use.
suppressPackageStartupMessages(library(landisR.fire))
stopifnot('output_dir' %in% names(formals(basefire_step)))
stopifnot('basefire_write_outputs' %in% getNamespaceExports('landisR.fire'))

out <- file.path(tempdir(), paste0('271_writer_test_',Sys.getpid()))
if (dir.exists(out)) unlink(out,recursive=TRUE)
land <- data.frame(cell_id=1:6,row=c(1,1,1,2,2,2),
                   col=c(1,2,3,1,2,3),active=c(TRUE,FALSE,TRUE,TRUE,TRUE,TRUE))
active <- land$cell_id[land$active]
severity <- data.frame(cell_id=active,severity=c(3,0,0,2,0),
                       disturbed=c(TRUE,FALSE,FALSE,TRUE,FALSE),
                       event_owner=c(1,0,0,2,0))
state <- list(cells=data.frame(cell_id=active,
                               last_fire_year=c(10,NA,6,10,NA),
                               last_severity=c(3,NA,1,2,NA)))
class(state) <- 'landisR_fire_state'
fixture <- list(
  simulation=list(year=10),
  fire_result=list(
    time=10L,
    events=list(list(metadata=list(event_id=1L),
                     severity=data.frame(cell_id=1L,severity=3)),
                list(metadata=list(event_id=2L),
                     severity=data.frame(cell_id=5L,severity=2))),
    raw_events=list(list(killed_cohort_rows=c(1L,2L)),
                    list(killed_cohort_rows=integer())),
    event_log=data.frame(event_id=1:2,time=10L,row=c(1L,2L),
                         column=c(1L,2L),sites_checked=c(2L,3L),
                         cohorts_killed=c(2L,0L),severity=c(3,2),
                         burned_sites=c(1L,1L)),
    summary=list(time=10L,number_fires=2L,total_sites_burned=2L,
                 region_sites=c(North=1L,South=1L)),
    fire_state=state,
    severity=severity
  ),
  prepared_profile=list(landscape=land,
                        model=list(prepared=list(size=list(cell_area=1))))
)
set.seed(471)
seed0 <- get('.Random.seed',envir=.GlobalEnv)
m <- basefire_write_outputs(fixture,out,formats=c('csv','rds'))
stopifnot(identical(seed0,get('.Random.seed',envir=.GlobalEnv)))
stopifnot(nrow(m)==7L,all(file.exists(m$path)))
stopifnot(file.exists(file.path(out,'year_0010','manifest.csv')))
ev <- utils::read.csv(file.path(out,'year_0010','fire_events.csv'))
rg <- utils::read.csv(file.path(out,'year_0010','fire_regions.csv'))
ce <- utils::read.csv(file.path(out,'year_0010','fire_event_cells.csv'))
stopifnot(nrow(ev)==2L,nrow(ce)==2L,nrow(rg)==2L,
          identical(sort(ce$cell_id),c(1L,5L)),
          identical(readRDS(file.path(out,'year_0010','fire_raw_events.rds')),
                    fixture$fire_result$raw_events))
err <- tryCatch({basefire_write_outputs(fixture,out);NULL},error=function(e)e)
stopifnot(inherits(err,'error'))
m2 <- basefire_write_outputs(fixture,out,formats=c('csv','rds'),overwrite=TRUE)
stopifnot(nrow(m2)==7L,identical(seed0,get('.Random.seed',envir=.GlobalEnv)))

# No fire: empty event CSV and valid history still written.
empty <- fixture
empty$fire_result$time <- 20L
empty$fire_result$summary$time <- 20L
empty$fire_result$summary$number_fires <- 0L
empty$fire_result$summary$total_sites_burned <- 0L
empty$fire_result$summary$region_sites[] <- 0L
empty$fire_result$events <- list()
empty$fire_result$raw_events <- list()
empty$fire_result$event_log <- empty$fire_result$event_log[FALSE,,drop=FALSE]
empty$fire_result$severity$severity[] <- 0
empty$fire_result$severity$disturbed[] <- FALSE
empty$fire_result$severity$event_owner[] <- 0
basefire_write_outputs(empty,out,formats='csv')
stopifnot(nrow(utils::read.csv(file.path(out,'year_0020','fire_events.csv')))==0L)

# GeoTIFF: verify row-major values, CRS and inactive mask when terra is available.
if (requireNamespace('terra',quietly=TRUE)) {
  ref <- terra::rast(nrows=2,ncols=3,xmin=100,xmax=130,ymin=50,ymax=70,crs='EPSG:4326')
  tifdir <- file.path(tempdir(),paste0('271_tif_test_',Sys.getpid()))
  if (dir.exists(tifdir)) unlink(tifdir,recursive=TRUE)
  x <- basefire_write_outputs(fixture,tifdir,formats='tif',raster_template=ref)
  stopifnot(nrow(x)==3L)
  r <- terra::rast(file.path(tifdir,'year_0010','rasters','fire_severity.tif'))
  z <- terra::values(r,mat=FALSE)
  stopifnot(length(z)==6L,isTRUE(all.equal(z[c(1,3,4,5,6)],c(3,0,0,2,0))),
            is.na(z[2]),terra::compareGeom(ref,r,stopOnError=FALSE))
  h <- terra::rast(file.path(tifdir,'year_0010','rasters','time_since_fire.tif'))
  stopifnot(isTRUE(all.equal(terra::values(h,mat=FALSE)[c(1,4,5)],c(0,4,0))))
  unlink(tifdir,recursive=TRUE)
}
cat('[PASS 271] CSV/RDS, empty-year, manifest, overwrite, RNG, raster (if terra)\n')
