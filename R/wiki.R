

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
  js$organisms
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
    return(NULL)
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
#' @param gene_dbs Names of databases to use
#'
#' @return A tibble with \code{term_id} and \code{gene_symbol}
fetch_wiki_pathway_genes_api <- function(pathways, gene_dbs = c("Ensembl", "HGNC", "HGNC Accession number", "Uniprot-TrEMBL")) {
  Database <- term_id <- TextLabel <- NULL

  pb <- progress::progress_bar$new(total = length(pathways))
  purrr::map_dfr(pathways, function(pathway) {
    pb$tick()
    js <- get_wiki_query("getPathway", stringr::str_glue("pwId={pathway}"))
    parse_wiki_gpml(js$pathway$gpml) |>
      tibble::add_column(term_id = pathway)
  }) |>
    dplyr::filter(Database %in% gene_dbs) |>
    dplyr::select(term_id, gene_symbol = TextLabel) |>
    dplyr::distinct()
}


#' Get functional term data from WikiPathways
#'
#' Download term information (pathway ID and name) and gene-pathway mapping
#' (gene symbol and pathway ID) from WikiPathways.
#'
#' @param species WikiPathways species designation, for example "Homo
#'   sapiens" for human. Full list of available species can be found using
#'   \code{fetch_wiki_species()}.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @import assertthat
#'
#' @examples
#' wiki_data <- fetch_wiki("Saccharomyces cerevisiae")
fetch_wiki <- function(species) {
  assert_that(is.string(species))

  terms <- fetch_wiki_pathways(species)
  assert_that(!is.null(terms), msg = stringr::str_glue("Error while extracting pathways from WikiPathways. Make sure species '{species}' is correct, use fetch_wiki_species()."))

  mapping <- fetch_wiki_pathway_genes_api(terms$term_id)
  list(
    terms = terms,
    mapping = mapping
  )
}
