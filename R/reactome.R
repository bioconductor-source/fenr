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
#' @param spec Reactome species
#'
#' @return A tibble with columns \code{term_id} and \code{term_name}
fetch_reactome_pathways <- function(spec) {
  # Binding variables from non-standard evaluation locally
  species <- NULL

  url <- "https://reactome.org/download/current/ReactomePathways.txt"
  assert_url_path(url)

  colms <- c("term_id", "term_name", "species")
  readr::read_tsv(url, col_names = colms, show_col_types = FALSE) |>
    dplyr::filter(species == spec)
}

#' Download term Ensembl gene ID mapping from Reactome
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


#' Get functional term data from Reactome
#'
#' Download term information (pathway ID and name) and gene-pathway mapping
#' (Ensembl gene ID and pathway ID) from Reactome.
#'
#' @param species Reactome species designation, for example "Homo
#'   sapiens" for human. Full list of available species can be found using
#'   \code{fetch_reactome_species()}.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#'
#' @examples
#' reactome_data <- fetch_reactome("Saccharomyces cerevisiae")
fetch_reactome <- function(species) {
  assert_that(!missing(species), msg = "Argument 'species' is missing.")
  assert_species(species, fetch_reactome_species)

  terms <- fetch_reactome_pathways(species)
  mapping <- fetch_reactome_ensembl_genes(species)
  list(
    terms = terms,
    mapping = mapping
  )
}
