#' A helper function to retrieve API data from Reactome
#'
#' This is used to quickly retrieve Reactome species. Alas, Reactome does not
#' provide API to download all pathways, so we need to download them directly
#' from a URL.
#'
#' @param path Path in the query
#'
#' @return A tibble with result
fetch_reactome_data <- function(path) {
  base_url <-  "https://reactome.org/ContentService"
  url <- file.path(base_url, path)
  assert_url_path(url)
  res <- httr::GET(url)
  jsonlite::fromJSON(rawToChar(res$content)) |>
    tibble::as_tibble()
}


#' List of available Reactome species
#'
#' @return A tibble with species names used by Reactome.
#' @export
#' @examples
#' re <- fetch_reactome_species()
fetch_reactome_species <- function() {
  # Binding variables from non-standard evaluation locally
  dbId <- displayName <- taxId <- NULL

  fetch_reactome_data("data/species/main") |>
    dplyr::select(db_id = dbId, designation = displayName, tax_id = taxId)
}


#' Download pathway data from Reactome
#'
#' @param tax_id Taxonomy ID of the species, a string
#'
#' @return A tibble with columns \code{term_id} and \code{term_name}
fetch_reactome_pathways <- function(tax_id) {
  # Binding variables from non-standard evaluation locally
  stId <- displayName <- NULL

  qry <- stringr::str_glue("data/schema/Pathway/min?species={tax_id}&page=1&offset=20000")
  fetch_reactome_data(qry) |>
    dplyr::select(term_id = stId, term_name = displayName)
}


#' Download term - Ensembl gene ID mapping from Reactome
#'
#' @details This function downloads one large file containing a mapping between
#'   Enxembl gene IDs and Reactome terms. This is significantly faster than
#'   \code{fetch_reactome_genes}, but if gene symbols are required, needs
#'   additional ID conversion.
#'
#' @param spec Reactome species
#'
#' @return A tibble with columns \code{gene_id} and \code{term_id}
fetch_reactome_ensembl_genes <- function(spec) {
  # Binding variables from non-standard evaluation locally
  species <- gene_id <- term_id <- NULL

  url <- "https://reactome.org/download/current/Ensembl2Reactome.txt"
  assert_url_path(url)
  colms <- c("gene_id", "term_id", "url", "event", "evidence", "species")
  readr::read_tsv(url, col_names = colms, show_col_types = FALSE) |>
    dplyr::filter(species == spec) |>
    dplyr::select(gene_id, term_id) |>
    dplyr::distinct()
}


#' Download term - gene symbol mapping from Reactome
#'
#' @details This function interrogates Reactome API to get term-gene mapping for
#'   all pathways. This is considerable slower than
#'   \code{fetch_reactome_ensembl_genes}. Warning, occasionally, for some
#'   pathways, Reactome does not return gene symbol - only UniProt accession
#'   number is available.
#'
#' @param pathways A character vector with Reactome patway IDs to get
#'   corresponding genes from.
#'
#' @return A tibble with columns\code{term_id}, \code{accession_number} and
#'   \code{gene_symbol}.
fetch_reactome_genes <- function(pathways) {
  identifier <- geneName <- gene_symbol <- databaseName <- NULL

  pb <- progress::progress_bar$new(total = length(pathways))
  purrr::map_dfr(pathways, function(pathway) {
    pb$tick()
    qry <- stringr::str_glue("data/participants/{pathway}/referenceEntities")
    fetch_reactome_data(qry) |>
      dplyr::filter(databaseName == "UniProt") |>
      dplyr::select(gene_symbol = geneName, accession_number = identifier) |>
      tidyr::unnest(gene_symbol) |>
      tibble::add_column(term_id = pathway, .before = 1)
  })
}


#' Get functional term data from Reactome
#'
#' Download term information (pathway ID and name) and gene-pathway mapping
#' (Ensembl gene ID or gene symbol and pathway ID) from Reactome.
#'
#' @details Reactome makes mapping between Ensembl ID and pathway ID available
#'   in form of one downloadable file. If `method = "file"` is set, this file
#'   will be downloaded and parsed. If `method = "api"` is set, then Reactome
#'   APIs will be interrogated for each pathway available. The first method is
#'   considerable faster and it returns Ensembl IDs. If gene symbols or UniProt
#'   accession numbers are needed, this will require additional conversion by
#'   the user. The second method is slower and it returns gene symbols.
#'
#' @param species Reactome species designation, for example "Homo sapiens" for
#'   human. Full list of available species can be found using
#'   \code{fetch_reactome_species()}.
#' @param method How to download the mapping. If 'file', then one mapping file
#'   provided by Reactome will be downloaded, if 'api', then Reactome API will
#'   be used. See details.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#'
#' @examples
#' reactome_data <- fetch_reactome("Saccharomyces cerevisiae")
fetch_reactome <- function(species, method = c("file", "api")) {
  method <- match.arg(method)
  assert_that(!missing(species), msg = "Argument 'species' is missing.")

  tax_id <- match_species(species, fetch_reactome_species, "tax_id")
  terms <- fetch_reactome_pathways(tax_id)
  if(method == "file") {
    mapping <- fetch_reactome_ensembl_genes(species)
  } else {
    mapping <- fetch_reactome_genes(terms$term_id)
  }

  list(
    terms = terms,
    mapping = mapping
  )
}
