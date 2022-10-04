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


test_that("Parsing KEGG flat file", {
  flat <- readr::read_file("../test_data/kegg_test.txt")
  expected <- readr::read_rds("../test_data/kegg_parseing_result.rds")
  parsed <- parse_kegg_genes(flat)
  expect_equal(parsed, expected)
})


test_that("Expected return from fetch_kegg_species", {
  expected_selection <- c("hsa", "mmu", "rno", "sce", "dme", "cel")
  spec <- fetch_kegg_species()
  expect_is(spec, "tbl")
  expect_true(all(expected_selection %in% spec$designation))
})


test_that("KEGG yeast mapping makes sense", {
  expected <- tibble::tribble(
    ~term_id, ~gene_symbol,
    "sce00010", "FBA1",
    "sce00010", "TDH3",
    "sce03030", "POL1",
    "sce04120", "UBI4"
  )
  mapping <- fetch_kegg_mapping(unique(expected$term_id), 1)
  expect_is(mapping, "tbl")
  merged <- expected |>
    dplyr::left_join(mapping, by = c("term_id", "gene_symbol")) |>
    tidyr::drop_na()
  expect_equal(nrow(expected), nrow(merged))
})
