test_that("Bioplanet mapping makes sense", {
  expected <- tibble::tribble(
    ~term_id, ~gene_symbol,
    "bioplanet_1025", "CDK1",
    "bioplanet_120", "IL1A",
    "bioplanet_1755", "RPL6",
    "bioplanet_1121", "POLR1A"
  )
  bp <- fetch_bp()
  expect_is(bp, "list")
  expect_setequal(names(bp), c("terms", "mapping"))
  mapping <- bp$mapping

  merged <- expected |>
    dplyr::left_join(mapping, by = c("term_id", "gene_symbol")) |>
    tidyr::drop_na()
  expect_equal(nrow(expected), nrow(merged))
})
