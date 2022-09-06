test_that("Incorrect species in fetch_reactome", {
  expect_error(fetch_reactome())
  expect_error(fetch_reactome(1243))
  expect_error(fetch_reactome("not a species"))
  expect_error(fetch_reactome("Homo sapiens", "blah"))
})

