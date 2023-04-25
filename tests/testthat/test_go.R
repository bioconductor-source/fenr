test_that("Incorrect URL in fetch_go_species", {
  expect_error(
    fetch_go_species("http://current.geneontology.org/products/pages/this_does_not_exist.html")
  )
})

test_that("Incorrect aruments in fetch_go", {
  expect_error(fetch_go())
  expect_error(fetch_go(species = "sgd", mart = "mart"))
})



test_that("Incorrect species in fetch_go", {
  expect_error(fetch_go(species = 1243))
  expect_error(fetch_go(species = "not a species"))
})


test_that("Incorrect mart in fetch_go", {
  expect_error(fetch_go(mart = 1))
  expect_error(fetch_go(mart = "mart"))
})


test_that("Processing OBO file", {
  expected_ids <- c("GO:0000001", "GO:0000022", "GO:0000232", "GO:1905121")
  expected_names <- c("mitochondrion inheritance", "mitotic spindle elongation",
                      "obsolete nuclear interphase chromosome", "mitotic spindle elongation")

  trms <- fetch_go_terms("../test_data/go_obo_test.txt")
  expect_s3_class(trms, "tbl")
  expect_equal(trms$term_id, expected_ids)
  expect_equal(trms$term_name, expected_names)
})


test_that("Expected return from fetch_go_species", {
  expected_selection <- c("goa_human", "mgi", "rgd", "sgd", "fb", "wb")
  spec <- fetch_go_species()
  expect_is(spec, "tbl")
  expect_true(all(expected_selection %in% spec$designation))
})


test_that("GO yeast gene mapping makes sense", {
  expected <- tibble::tribble(
    ~term_id, ~gene_symbol,
    "GO:0000166", "POL1",
    "GO:0006096", "FBA1",
    "GO:0005199", "CWP2",
    "GO:0004365", "TDH3"
  )
  mapping <- fetch_go_genes_go("sgd") |>
    dplyr::select(term_id, gene_symbol) |>
    dplyr::distinct()
  expect_is(mapping, "tbl")
  merged <- expected |>
    dplyr::left_join(mapping, by = c("term_id", "gene_symbol")) |>
    tidyr::drop_na()
  expect_equal(nrow(expected), nrow(merged))
})
