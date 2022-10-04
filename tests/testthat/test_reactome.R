test_that("Incorrect species in fetch_reactome", {
  expect_error(fetch_reactome())
  expect_error(fetch_reactome(1243))
  expect_error(fetch_reactome("not a species"))
  expect_error(fetch_reactome("Homo sapiens", "blah"))
})


test_that("Expected return from fetch_reactome_species", {
  expected_selection <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus",
                          "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans")
  spec <- fetch_reactome_species()
  expect_is(spec, "tbl")
  expect_true(all(expected_selection %in% spec$designation))
})


test_that("Reactome yeast mapping makes sense", {
  species <- "Saccharomyces cerevisiae"
  tax_id <- "4932"

  expected <- tibble::tribble(
    ~term_id, ~gene_id,
    "R-SCE-70171", "YJL052W",
    "R-SCE-70171", "YGR192C",
    "R-SCE-68952", "YNL102W",
    "R-SCE-983168", "YKL022C"
  )

  mapping <- fetch_reactome_ensembl_genes(species)
  expect_is(mapping, "tbl")
  merged <- expected |>
    dplyr::left_join(mapping, by = c("term_id", "gene_id")) |>
    tidyr::drop_na()
  expect_equal(nrow(expected), nrow(merged))
})
