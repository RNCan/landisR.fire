test_that('fire_event stores cell-level severity', {
  e <- fire_event(
    year = 10, ignition_cell = 2, cells = c(2,3,4),
    severity = c(2,4,5), source_cell = c(NA,2,3)
  )
  expect_s3_class(e, 'landisR_fire_event')
  expect_equal(e$severity$severity, c(2,4,5))
  expect_equal(e$spread$source_cell, c(NA,2,3))
})

test_that('scalar severity expands across cells', {
  e <- fire_event(0, 1, c(1,2,3), severity = 3)
  expect_equal(e$severity$severity, c(3,3,3))
})

test_that('duplicate burned cells are rejected', {
  expect_error(fire_event(0, 1, c(1,1), severity = 2), 'unique')
})
