test_that("Incorrect species in fetch_kegg", {
  expect_error(fetch_kegg())
  expect_error(fetch_kegg(1243))
  expect_error(fetch_kegg("not a species"))
})

test_that("Incorrect batch size in fetch_kegg", {
  expect_error(fetch_kegg("sce", "aha!"))
  expect_error(fetch_kegg("sce", -1))
  expect_error(fetch_kegg("sce", 0))
  expect_error(fetch_kegg("sce", 1000))
})
