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


test_that("Correct response from fetch_kegg", {
  species <- "sce"

  expected_terms <- tibble::tribble(
    ~term_id, ~term_name,
    "sce01100",  "Metabolic pathways - Saccharomyces cerevisiae (budding yeast)",
    "sce01230", "Biosynthesis of amino acids - Saccharomyces cerevisiae (budding yeast)",
    "sce03010", "Ribosome - Saccharomyces cerevisiae (budding yeast)",
    "sce00052", "Galactose metabolism - Saccharomyces cerevisiae (budding yeast)"
  )

  expected_mapping <- tibble::tribble(
    ~term_id, ~gene_symbol,
    "sce00010", "FBA1",
    "sce00010", "TDH3",
    "sce03030", "POL1",
    "sce04120", "UBI4"
  )

  re <- fetch_kegg("sce")
  test_fetched_structure(re)
  test_terms(re$terms, expected_terms)
  test_mapping(re$mapping, expected_mapping, "gene_symbol")
})
