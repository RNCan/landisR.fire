test_that('empirical size distribution is reproducible and respects support', {
  d <- fire_size_empirical(c(10, 100, 1000), c(0.7, 0.2, 0.1))
  set.seed(4357)
  a <- d(50)
  set.seed(4357)
  b <- d(50)
  expect_identical(a, b)
  expect_true(all(a %in% c(10, 100, 1000)))
})

test_that('empirical probabilities are normalized', {
  d <- fire_size_empirical(c(1, 2), c(2, 1))
  expect_equal(attr(d, 'probabilities'), c(2/3, 1/3))
})
