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


test_that("Expected return from fetch_wiki_species", {
  expected_selection <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus",
                          "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans")
  spec <- fetch_wiki_species()
  expect_is(spec, "tbl")
  expect_true(all(expected_selection %in% spec$designation))
})


test_that("WikiPathways yeast mapping makes sense", {
  databases = c("Ensembl", "Entrez Gene", "HGNC", "HGNC Accession number", "Uniprot-TrEMBL")
  types = c("GeneProduct", "Protein", "Rna", "RNA")

  expected <- tibble::tribble(
    ~term_id, ~gene_symbol,
    "WP253", "FBA1",
    "WP253", "TDH3",
    "WP13", "POL1",
    "WP13", "UBI4"
  )
  mapping <- fetch_wiki_pathway_genes_api(expected$term_id, databases, types)
  expect_is(mapping, "tbl")
  merged <- expected |>
    dplyr::left_join(mapping, by = c("term_id", "gene_symbol" = "text_label")) |>
    tidyr::drop_na()
  expect_equal(nrow(expected), nrow(merged))
})
