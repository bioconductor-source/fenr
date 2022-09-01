test_that("Incorrect species in fetch_wiki", {
  expect_error(fetch_wiki())
  expect_error(fetch_wiki(1243))
  expect_error(fetch_wiki("not a species"))
})

test_that("Parsing GPLM", {
  expected <- tibble::tibble(
    TextLabel = c("STR3", "CYS4", "STR2"),
    Type = c("GeneProduct", "GeneProduct", "GeneProduct"),
    Database = c("SGD", "SGD", "SGD"),
    ID = c("S000003152", "S000003387", "S000003891")
  )

  gpml <- readLines("../test_data/wiki_test.gpml") |>
    paste(collapse = "\n")
  res <- parse_wiki_gpml(gpml)

  expect_equal(res, expected)
})
