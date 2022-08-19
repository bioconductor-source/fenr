#' Download pathway data from Reactome
#'
#' @param spec Reactome species
#'
#' @return A tibble with columns \code{term_id} and \code{term_name}
fetch_reactome_pathways <- function(spec) {
  # Binding variables from non-standard evaluation locally
  species <- spec <- NULL

  u <- "https://reactome.org/download/current/ReactomePathways.txt"
  stopifnot(url_exists(u))
  colms <- c("term_id", "term_name", "species")
  readr::read_tsv(u, col_names = colms, show_col_types = FALSE) |>
    dplyr::filter(species == spec)
}

#' Download term Ensembl gene ID mapping from Reactome
#'
#' @param spec Reactome species
#'
#' @return A tibble with columns \code{gene_id} and \code{term_id}
fetch_reactome_ensembl_genes <- function(spec) {
  # Binding variables from non-standard evaluation locally
  species <- spec <- gene_id <- term_id <- NULL

  u <- "https://reactome.org/download/current/Ensembl2Reactome.txt"
  stopifnot(url_exists(u))
  colms <- c("gene_id", "term_id", "url", "event", "evidence", "species")
  readr::read_tsv(u, col_names = colms, show_col_types = FALSE) |>
    dplyr::filter(species == spec) |>
    dplyr::select(gene_id, term_id) |>
    dplyr::distinct()
}

#' List of available Reactome species
#'
#' @return A character vector with species names used by Reactome.
#' @export
fetch_reactome_species <- function() {
  # Binding variables from non-standard evaluation locally
  species <- NULL

  u <- "https://reactome.org/download/current/Ensembl2Reactome.txt"
  stopifnot(url_exists(u))
  colms <- c("gene_id", "term_id", "url", "event", "evidence", "species")
  readr::read_tsv(u, col_names = colms, show_col_types = FALSE) |>
    dplyr::pull(species) |>
    unique()
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
#' @import assertthat
#'
#' @examples
#' \dontrun{
#' reactome_data <- fetch_reactome("Homo sapiens")
#' }
fetch_reactome <- function(species) {
  assert_that(is.string(species))

  terms <- fetch_reactome_pathways(species)
  assert_that(nrow(terms) > 0, msg = stringr::str_glue("No pathways found from Reactome. Make sure species '{species}' is correct, use fetch_reactome_species()."))

  mapping <- fetch_reactome_ensembl_genes(species)
  list(
    terms = terms,
    mapping = mapping
  )
}
