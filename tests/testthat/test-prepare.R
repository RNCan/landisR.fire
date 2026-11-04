test_that('fire_prepare calls component preparation once and stores caches', {
  landscape <- expand.grid(row = 1:3, col = 1:3)
  landscape$cell_id <- seq_len(nrow(landscape))
  landscape$active <- TRUE

  calls <- 0L

  spread <- fire_spread_custom(
    function(..., prepared = NULL) prepared,
    prepare = function(landscape, fire_state, topology, context, model) {
      calls <<- calls + 1L
      list(
        n = nrow(landscape),
        degree = rowSums(topology$neighbors > 0)
      )
    }
  )

  f <- function(...) NULL

  model <- fire_model(
    ignition = f,
    size = f,
    spread = spread,
    burn_probability = f,
    severity = f,
    mortality = f
  )

  prepared <- fire_prepare(model, landscape)

  expect_s3_class(prepared, 'landisR_fire_prepared_model')
  expect_equal(calls, 1L)
  expect_equal(prepared$prepared$spread$n, 9)
  expect_equal(length(prepared$prepared$spread$degree), 9)
  expect_s3_class(prepared$topology, 'landisR_fire_topology')
})

test_that('components without prepare functions have NULL caches', {
  landscape <- data.frame(
    cell_id = 1:4,
    row = c(1, 1, 2, 2),
    col = c(1, 2, 1, 2)
  )

  f <- function(...) NULL
  model <- fire_model(f, f, f, f, f, f)
  prepared <- fire_prepare(model, landscape)

  expect_true(all(vapply(prepared$prepared, is.null, logical(1))))
})
