# 216B-public-api-v1

test_that('basefire_read is a public parser entry point', {
  expect_error(
    basefire_read(tempfile('missing-basefire-')),
    'Base Fire parameter file not found'
  )
})

test_that('Base Fire public workflow is intentionally minimal', {
  exports <- getNamespaceExports('landisR.fire')

  basefire_exports <- sort(
    exports[grepl('^basefire_', exports)]
  )

  expect_identical(
    basefire_exports,
    sort(c("basefire_prepare", "basefire_read", "basefire_step", "basefire_write_outputs"))
  )
})
