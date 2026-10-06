test_that('fire_model stores replaceable components', {
  f <- function(...) NULL
  m <- fire_model(f, f, f, f, f, name = 'x')
  expect_s3_class(m, 'landisR_fire_model')
  expect_identical(m$name, 'x')
  expect_true(is.function(m$spread))
})
