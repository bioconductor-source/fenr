species <- "Saccharomyces cerevisiae"
tax_id <- "4932"

expected_pathways <- tibble::tribble(
  ~term_id, ~term_name,
  "R-SCE-6798695", "Neutrophil degranulation",
  "R-SCE-983168", "Antigen processing: Ubiquitination & Proteasome degradation",
  "R-SCE-204005", "COPII-mediated vesicle transport",
  "R-SCE-5668541", "TNFR2 non-canonical NF-kB pathway"
)

expected_ensembl <- tibble::tribble(
  ~term_id, ~gene_id,
  "R-SCE-70171", "YJL052W",
  "R-SCE-70171", "YGR192C",
  "R-SCE-68952", "YNL102W",
  "R-SCE-983168", "YKL022C"
)

expected_gene_association <- tibble::tribble(
  ~term_id, ~gene_symbol,
  "R-SCE-6793739", "FBRL",
  "R-SCE-9749345", "CDT1",
  "R-SCE-1252249", "ATPB",
  "R-SCE-9749381", "CDC6",
)

expected_api <- tibble::tribble(
  ~term_id, ~gene_symbol,
  "R-SCE-68952", "POL12",
  "R-SCE-983168", "CDC16",
  "R-HSA-114608", "CFD",
  "R-MMU-352230", "Slc7a5"
)

##################################################

test_that("Incorrect species in fetch_reactome", {
  expect_error(fetch_reactome())
  expect_error(fetch_reactome(1243))
  expect_error(fetch_reactome("not a species"))
  expect_error(fetch_reactome("Homo sapiens", "blah"))
})

test_that("Expected behaviour from a non-responsive server", {
  url <- get_reactome_url()
  options(REACTOME_BASE_URL = SERVER_500)
  test_unresponsive_server(fetch_reactome_species)
  test_unresponsive_server(fetch_reactome_pathways, tax_id = tax_id)
  test_unresponsive_server(fetch_reactome, species = species)
  options(REACTOME_BASE_URL = url)

  url <- get_reactome_ensembl_file()
  options(REACTOME_ENSEMBL_FILE = SERVER_500)
  test_unresponsive_server(fetch_reactome_ensembl_genes, spec = species, use_cache = FALSE)
  options(REACTOME_ENSEMBL_FILE = url)

  url <- get_reactome_gaf_file()
  options(REACTOME_GAF_FILE = SERVER_500)
  test_unresponsive_server(fetch_reactome_gene_association, tax_id = tax_id, use_cache = FALSE)
  options(REACTOME_GAF_FILE = url)
})

test_that("Expected return from fetch_reactome_species", {
  expected_selection <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus",
                          "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans")
  spec <- fetch_reactome_species(on_error = "warn")
  if(!is.null(spec)) {
    expect_is(spec, "tbl")
    expect_true(all(expected_selection %in% spec$designation))
  }
})


test_that("Correct Reactome Ensembl from yeast", {
  re <- fetch_reactome(species, source = "ensembl", use_cache = TRUE, on_error = "warn")
  if(!is.null(re)) {
    expect_is(re, "list")
    expect_length(re, 2)
    expect_named(re)
    expect_equal(names(re), c("terms", "mapping"))

    # Check pathways
    paths <- re$terms
    expect_is(paths, "tbl")

    merged <- expected_pathways |>
      dplyr::left_join(paths, by = c("term_id", "term_name")) |>
      tidyr::drop_na()
    expect_equal(nrow(expected_pathways), nrow(merged))

    # Check mapping
    mapping <- re$mapping
    expect_is(mapping, "tbl")
    merged <- expected_ensembl |>
      dplyr::left_join(mapping, by = c("term_id", "gene_id")) |>
      tidyr::drop_na()
    expect_equal(nrow(expected_ensembl), nrow(merged))
  }
})


test_that("Correct Reactome gene association from yeast", {
  re <- fetch_reactome(species, source = "gene_association", use_cache = TRUE, on_error = "warn")
  if(!is.null(re)) {
    expect_is(re, "list")
    expect_length(re, 2)
    expect_named(re)
    expect_equal(names(re), c("terms", "mapping"))

    # Check mapping
    mapping <- re$mapping
    expect_is(mapping, "tbl")
    merged <- expected_gene_association |>
      dplyr::left_join(mapping, by = c("term_id", "gene_symbol")) |>
      tidyr::drop_na()
    expect_equal(nrow(expected_gene_association), nrow(merged))
  }
})


test_that("Correct Reactome API from yeast", {
  # Doing this on all pathways takes forever
  mapping <- fetch_reactome_api_genes(expected_api$term_id, on_error = "warn")
  if(!is.null(mapping)) {
    expect_is(mapping, "tbl")
    merged <- expected_api |>
      dplyr::left_join(mapping, by = c("term_id", "gene_symbol")) |>
      tidyr::drop_na()
    expect_equal(nrow(expected_api), nrow(merged))
  }
})
