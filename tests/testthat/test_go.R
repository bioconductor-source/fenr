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


test_that("Processing OBO file", {
  expected_ids <- c("GO:0000001", "GO:0000022", "GO:0000232")
  expected_names <- c("mitochondrion inheritance", "mitotic spindle elongation",
                      "obsolete nuclear interphase chromosome")
  obo <- readr::read_lines("../test_data/go_obo_test.txt")
  parsed <- parse_obo_file(obo)
  ids <- parsed |>
    dplyr::filter(key == "id") |>
    dplyr::pull(value)
  names <- parsed |>
    dplyr::filter(key == "name") |>
    dplyr::pull(value)
  expect_vector(ids)
  expect_vector(names)
  expect_equal(ids, expected_ids)
  expect_equal(names, expected_names)
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
  mapping <- fetch_go_genes_go("sgd")
  expect_is(mapping, "tbl")
  merged <- expected |>
    dplyr::left_join(mapping, by = c("term_id", "gene_symbol")) |>
    tidyr::drop_na()
  expect_equal(nrow(expected), nrow(merged))
})
