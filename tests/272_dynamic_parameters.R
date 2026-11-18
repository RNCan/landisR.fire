# Package self-tests. No LANDIS-II project or full simulation required.
suppressPackageStartupMessages(library(landisR.fire))
stopifnot(all(c('basefire_parameter_schedule', 'basefire_update_parameters',
                'basefire_step_dynamic') %in% getNamespaceExports('landisR.fire')))
private <- function(nm) get(nm, envir=asNamespace('landisR.fire'))
param_fun <- private('.basefire_parameters')
profile_fun <- private('.basefire_profile')
prepare_fun <- private('.basefire_prepare_profile')

regions <- data.frame(name=c('North','South'), map_code=c(1L,2L),
  mean_size=c(100,250), min_size=c(10,20), max_size=c(700,900),
  ignition_probability=c(0.01,0.02), fire_spread_age=c(50L,120L))
params <- param_fun(
  timestep=10L,
  species=data.frame(species_code='BlackSpruce', fire_tolerance=2L),
  regions=regions,
  initial_fire_regions='initial_regions.tif',
  dynamic_regions=data.frame(year=20L, map_name='regions_future.tif'),
  damage=data.frame(max_age=c(0.4,1.0),
                    severity_tolerance_diff=c(-2L,2L)),
  map_names_template='fire-{timestep}.tif',
  log_file='fire.log', summary_log_file='fire_summary.csv'
)
land <- data.frame(cell_id=1:4, row=c(1L,1L,2L,2L),
  col=c(1L,2L,1L,2L), active=rep(TRUE,4),
  fire_region_map_code=c(1L,1L,2L,2L),
  max_cohort_age=c(80,70,60,50))
map_reader <- function(map_name, current_time, landscape) {
  landscape$fire_region_map_code
}
ctx <- list(start_time=0L, cell_area=1,
            species_longevity=c(BlackSpruce=180),
            dynamic_region_map_reader=map_reader)
p <- prepare_fun(profile_fun(params),landscape=land,
                 fire_state=NULL,context=ctx,topology=NULL)

changes <- data.frame(year=c(10L,20L,20L,30L),
                      region_name=c('North','North','South','North'),
                      ignition_probability=c(0.03, NA, 0.05, 0.04),
                      mean_size=c(NA, 140, NA, NA),
                      fire_spread_age=c(NA,40L,150L,NA))
schedule <- basefire_parameter_schedule(changes)
set.seed(272L)
seed0 <- .Random.seed
p0 <- basefire_update_parameters(p,schedule,current_time=0L)
stopifnot(identical(seed0,.Random.seed),
          identical(p0$fire_state,p$fire_state),
          identical(p0$parameters$regions,p$parameters$regions))
p10 <- basefire_update_parameters(p0,schedule,current_time=10L)
stopifnot(identical(seed0,.Random.seed),
          p10$parameters$regions$ignition_probability[1] == 0.03,
          p10$parameter_cache$ignition_probability[1] == 0.03,
          identical(p10$fire_state,p$fire_state),
          length(p10$context$parameter_change_log$parameter)==1L)
p15 <- basefire_update_parameters(p10,schedule,current_time=15L)
stopifnot(identical(p15$parameter_cache,p10$parameter_cache))
p20 <- basefire_update_parameters(p15,schedule,current_time=20L)
stopifnot(p20$parameters$regions$mean_size[1] == 140,
          p20$parameters$regions$fire_spread_age[1] == 40,
          p20$parameters$regions$ignition_probability[2] == 0.05,
          p20$model$prepared$size$mean_size[1] == 140,
          identical(p20$fire_state,p$fire_state),
          identical(p20$context$dynamic_region_map_reader, map_reader),
          identical(p20$landscape$fire_region_map_code,
                    p$landscape$fire_region_map_code),
          identical(p20$parameter_cache$dynamic_year, 20L),
          identical(seed0,.Random.seed))
p20again <- basefire_update_parameters(p20,schedule,current_time=20L)
stopifnot(identical(p20again$context$parameter_change_log,
                    p20$context$parameter_change_log))
p30 <- basefire_update_parameters(p20again,schedule,current_time=30L)
stopifnot(p30$parameters$regions$ignition_probability[1] == 0.04)

geometry <- private(".basefire_shared_geometry")(land)
shared <- prepare_fun(profile_fun(params), landscape = land, context = ctx,
                      topology = geometry$active, site_topology = geometry$sites)
shared20 <- basefire_update_parameters(shared, schedule, current_time = 20L)
stopifnot(identical(shared20$model$prepared, p20$model$prepared),
          identical(shared20$topology, geometry$active),
          identical(seed0, .Random.seed))

# Map-code keyed updates are supported, but region identifiers remain fixed.
bycode <- basefire_parameter_schedule(data.frame(
  year=10L, map_code=2L, ignition_probability=0.11))
px <- basefire_update_parameters(p,bycode,current_time=10L)
stopifnot(px$parameters$regions$ignition_probability[2] == 0.11)

# Guard rails: invalid values, duplicate keys, unknown region, backwards time,
# and attempts to swap a schedule mid-run must all fail before changing state.
fail <- function(expr) inherits(tryCatch(force(expr),error=identity),'error')
stopifnot(fail(basefire_parameter_schedule(data.frame(
  year=10,region_name='North',ignition_probability=1.3))),
  fail(basefire_parameter_schedule(data.frame(
  year=c(10,10),region_name=c('North','North'),mean_size=c(1,2)))),
  fail(basefire_update_parameters(p,data.frame(
  year=10,region_name='West',mean_size=10),10L)),
  fail(basefire_update_parameters(p20,schedule,10L)),
  fail(basefire_update_parameters(p20,bycode,30L)),
  identical(seed0,.Random.seed))
cat('[PASS 272] preparation + cache refresh + time + idempotence + validation + RNG\n')
