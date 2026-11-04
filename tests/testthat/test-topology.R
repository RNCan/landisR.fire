test_that('regular-grid topology precomputes first and second directional rings', {
  landscape <- expand.grid(row = 1:3, col = 1:3)
  landscape$cell_id <- seq_len(nrow(landscape))
  landscape$active <- TRUE

  topo <- fire_grid_topology(landscape)

  expect_s3_class(topo, 'landisR_fire_topology')
  expect_equal(dim(topo$neighbors), c(9, 8))
  expect_equal(dim(topo$second_neighbors), c(9, 8))

  center <- which(topo$row == 2 & topo$col == 2)
  expect_equal(sum(topo$neighbors[center, ] > 0), 8)
  expect_equal(sum(topo$second_neighbors[center, ] > 0), 0)

  corner <- which(topo$row == 1 & topo$col == 1)
  expect_equal(sum(topo$neighbors[corner, ] > 0), 3)
})

test_that('inactive cells are absent from topology and behave as gaps', {
  landscape <- expand.grid(row = 1:3, col = 1:3)
  landscape$cell_id <- seq_len(nrow(landscape))
  landscape$active <- TRUE
  landscape$active[landscape$row == 2 & landscape$col == 2] <- FALSE

  topo <- fire_grid_topology(landscape)

  expect_equal(length(topo$cell_id), 8)
  expect_false(5 %in% topo$cell_id)
})
