#' REST API query for WikiPathways
#'
#' @param query Query, e.g. `listPathways`
#' @param parameters Parameters, e.g. `organism=Homo sapiens`
#'
#' @return JSON with response
get_wiki_query <- function(query, parameters = NULL) {
  pars <- paste(c(parameters, "format=json"), collapse = "&")
  qry <- stringr::str_glue("https://webservice.wikipathways.org/{query}?{pars}") |>
    stringr::str_replace_all("\\s", "%20")
  assert_http_file(qry)
  resp <- httr::GET(url = qry)
  jsonlite::fromJSON(rawToChar(resp$content))
}


#' List of available WikiPathways species
#'
#' @return A character vector with species names used by WikiPathways.
#' @export
#'
#' @examples
#' spec <- fetch_wiki_species()
fetch_wiki_species <- function() {
  js <- get_wiki_query("listOrganisms")
  tibble::tibble(designation = js$organisms)
}

#' Download pathway data from WikiPathways
#'
#' @param species Species name recognised by WikiPathways (see \code{fetch_wiki_species()})
#'
#' @return A tibble with columns \code{gene_id} and \code{term_id}
fetch_wiki_pathways <- function(species) {
  id <- name <- NULL

  js <- get_wiki_query("listPathways", stringr::str_glue("organism={species}"))
  if(js[1] == "error")
    stop(stringr::str_glue("Cannot retrieve pathways from WikiPathways for species {species}."))
  js$pathways |>
    tibble::as_tibble() |>
    dplyr::select(term_id = id, term_name = name)
}


#' Helper function for `parse_wiki_gpml`
#'
#' @param s Input string
#' @param key Key to extract
#'
#' @return Value of the key
extract_key <- function(s, key) {
  stringr::str_extract(s, paste0('(?<=', key, '=")(.+?)(?=")'))
}

#' Parse GPML string provided by WikiPathways
#'
#' @param gpml GMPL string returned by WikiPathways query.
#' @param keys Keys to be extracted.
#'
#' @return A tibble with selected keys and values
parse_wiki_gpml <- function(gpml, keys = c("TextLabel", "Type", "Database", "ID")) {
  gpml |>
    stringr::str_replace_all("\\n", "") |>
    stringr::str_extract_all("<DataNode.+?</DataNode>") |>
    purrr::pluck(1) |>
    purrr::map_dfr(function(s) {
      purrr::map_chr(keys, ~extract_key(s, .x)) |>
      purrr::set_names(keys)
    })
}


#' Download WikiPathways gene-pathway mapping using API
#'
#' @param pathways A character vector with pathway names.
#' @param databases Names of databases to use
#' @param types Names of types to use
#'
#' @return A tibble with \code{term_id} and \code{gene_symbol}
fetch_wiki_pathway_genes_api <- function(pathways, databases = NULL, types = NULL) {
  term_id <- TextLabel <-ID <- Type <- Database <- database <- type <- NULL

  pb <- progress::progress_bar$new(total = length(pathways))
  res <- purrr::map_dfr(pathways, function(pathway) {
    pb$tick()
    js <- get_wiki_query("getPathway", stringr::str_glue("pwId={pathway}"))
    parse_wiki_gpml(js$pathway$gpml) |>
      tibble::add_column(term_id = pathway)
  }) |>
    dplyr::select(term_id, text_label = TextLabel, id = ID, type = Type, database = Database) |>
    dplyr::distinct()
  if(!is.null(databases)) {
    res <- res |>
      dplyr::filter(database %in% databases)
  }
  if(!is.null(types)) {
    res <- res |>
      dplyr::filter(type %in% types)
  }

  res
}


#' Get functional term data from WikiPathways
#'
#' Download term information (pathway ID and name) and gene-pathway mapping
#' (gene symbol and pathway ID) from WikiPathways.
#'
#' @details WikiPathways contain mapping between pathways and a variety of
#'   entities from various databases. Typically a gene symbol is returned in
#'   column \code{text_label} and some sort of ID in column \code{id}, but this
#'   depends on the species and databases used. For gene/protein enrichment,
#'   these should be filtered to contain gene symbols only. This can be done by
#'   selecting a desired databases and types. The default values for parameters
#'   \code{databases} and \code{types} attempt to select inforation from generic
#'   databases, but there are organism-specific databases not included in the
#'   selection. We suggest to run this function with \code{databases = NULL,
#'   types = NULL} to see what types and databases are available before making
#'   selection.
#'
#' @param species WikiPathways species designation, for example "Homo sapiens"
#'   for human. Full list of available species can be found using
#'   \code{fetch_wiki_species()}.
#' @param databases A character vector with a list of databases to pre-filter mapping
#'   data. See details. Full result will be returned if NULL.
#' @param types A character vector with a list of types to pre-filter mapping data.
#'   See details. Full result will be returned if NULL.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @import assertthat
#'
#' @examples
#' wiki_data <- fetch_wiki("Saccharomyces cerevisiae")
fetch_wiki <- function(species,
    databases = c("Ensembl", "Entrez Gene", "HGNC", "HGNC Accession number", "Uniprot-TrEMBL"),
    types = c("GeneProduct", "Protein", "Rna", "RNA")
                       ) {
  assert_that(!missing(species), msg = "Argument 'species' is missing.")
  assert_species(species, fetch_wiki_species)

  terms <- fetch_wiki_pathways(species)
  mapping <- fetch_wiki_pathway_genes_api(terms$term_id, databases, types)
  list(
    terms = terms,
    mapping = mapping
  )
}
