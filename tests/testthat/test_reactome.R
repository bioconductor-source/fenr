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


test_that("Reactome Ensembl yeast mapping makes sense", {
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


test_that("Reactome API mapping makes sense", {
  expected <- tibble::tribble(
    ~term_id, ~accession_number, ~gene_symbol,
    "R-SCE-68952", "P38121", "POL12",
    "R-SCE-983168", "P09798", "CDC16",
    "R-HSA-114608", "P20160", "AZU1",
    "R-MMU-352230", "Q9Z127", "Slc7a5"
  )

  mapping <- fetch_reactome_genes(expected$term_id)
  expect_is(mapping, "tbl")
  merged <- expected |>
    dplyr::left_join(mapping, by = c("term_id", "accession_number", "gene_symbol")) |>
    tidyr::drop_na()
  expect_equal(nrow(expected), nrow(merged))
})
