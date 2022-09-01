test_that("Incorrect URL in fetch_go_species", {
  expect_error(
    fetch_go_species("http://current.geneontology.org/products/pages/this_does_not_exist.html")
  )
})

test_that("Incorrect species in fetch_go_from_go", {
  expect_error(fetch_go_from_go())
  expect_error(fetch_go_from_go(1243))
  expect_error(fetch_go_from_go("not a species"))
})

test_that("Incorrect mart in fetch_go_from_bm", {
  expect_error(fetch_go_from_bm())
  expect_error(fetch_go_from_bm(1))
  expect_error(fetch_go_from_bm("mart"))
})
