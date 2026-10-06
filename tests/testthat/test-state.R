test_that('fire_state stores independent fire history', {
  s <- fire_state(1:5)
  expect_s3_class(s, 'landisR_fire_state')
  expect_equal(nrow(s$cells), 5)
  expect_true(all(is.na(s$cells$last_fire_year)))
})
