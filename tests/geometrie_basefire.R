legacy <- getFromNamespace(".basefire_build_site_topology", "landisR.fire")
prepare <- getFromNamespace(".basefire_shared_geometry", "landisR.fire")
grid <- expand.grid(row = -2L:2L, col = 2L:7L)
grid <- grid[-c(3L, 12L), ]
set.seed(9)
grid <- grid[sample.int(nrow(grid)), ]
grid$cell_id <- seq_len(nrow(grid)) * 3L
grid$active <- seq_len(nrow(grid)) %% 3L != 0L
for (landscape in list(grid, transform(grid, active = TRUE), grid[1L, ])) {
  topology <- landisR.fire::fire_grid_topology(landscape)
  reference <- legacy(landscape, topology)
  actual <- prepare(landscape)
  stopifnot(identical(actual$sites, reference),
            identical(actual$active$cell_id, topology$cell_id),
            identical(actual$active$row, topology$row),
            identical(actual$active$col, topology$col))
}
invalid <- grid
invalid[2L, c("row", "col")] <- invalid[1L, c("row", "col")]
stopifnot(inherits(try(prepare(invalid), silent = TRUE), "try-error"))
