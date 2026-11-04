test_that('fire_model stores six replaceable components', {
  f <- function(...) NULL
  m <- fire_model(f, f, f, f, f, f, name = 'x')

  expect_s3_class(m, 'landisR_fire_model')
  expect_identical(m$name, 'x')
  expect_true(is.function(m$spread))
  expect_true(is.function(m$burn_probability))
  expect_identical(attr(m$spread, 'fire_component_class'), 'spread')
  expect_identical(
    attr(m$burn_probability, 'fire_component_class'),
    'burn_probability'
  )
})

test_that('component type mismatches are rejected', {
  f <- function(...) NULL
  spread <- fire_spread_custom(f)
  bad_burn <- spread

  expect_error(
    fire_model(f, f, spread, bad_burn, f, f),
    'incompatible'
  )
})
