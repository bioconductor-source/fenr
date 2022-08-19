GAF_COLUMNS <- c(
  "db",
  "db_id",
  "symbol",
  "qualifier",
  "go_term",
  "db_ref",
  "evidence",
  "from",
  "aspect",
  "db_object_name",
  "db_object_synonym",
  "db_object_type",
  "taxon",
  "date",
  "assigned_by",
  "annotation_extension",
  "form_id"
)


#' Download GO term descriptions
#'
#' @param obo_file A URL or a local file name containing GO ontology, in OBO format.
#'
#' @return A tibble with term_id and term_name.
fetch_go_terms <- function(obo_file = "http://purl.obolibrary.org/obo/go.obo") {
  # Binding variables from non-standard evaluation locally
  term_id <- term_name <- NULL

  stopifnot(url_exists(obo_file))
  go <- ontologyIndex::get_ontology(obo_file, extract_tags = "minimal")
  tibble::tibble(
    term_id = go$id,
    term_name = go$name
  )
}

#' Find all species available from geneontology.org
#'
#' This function attempts to scrape HTML web page containing a table of available species and corresponding file names. If the structure of the page changes one day and the function stops working, go to \url{http://current.geneontology.org/products/pages/downloads.html} and check file names. The species designation used in this package is the GAF file name without extension (e.g. for a file \file{goa_chicken.gaf} the designation is \file{goa_chicken}).
#'
#' @return A tibble with columns \code{species} and \code{designation}.
#' @export
fetch_go_species <- function() {
  # Binding variables from non-standard evaluation locally
  species <- designation <- `Species/Database` <- File <- NULL

  u <- "http://current.geneontology.org/products/pages/downloads.html"
  stopifnot(url_exists(u))
  u <- httr::GET(u) |>
    httr::content("text") |>
    XML::readHTMLTable(as.data.frame = TRUE)
  u[[1]] |>
    tibble::as_tibble() |>
    dplyr::mutate(
      species = `Species/Database` |>
        stringr::str_replace_all("\\n", "-") |>
        stringr::str_replace_all("\\s\\s+", " ") |>
        stringr::str_replace_all("(\\S)-", "\\1"),
      designation = File |>
        stringr::str_remove("\\..*$")
    ) |>
    dplyr::select(species, designation)
}


#' Download GO term gene mapping from geneontology.org
#'
#' @param species Species designation. Base file name for species file under
#'   \url{http://current.geneontology.org/annotations}. Examples are
#'   \file{goa_human} for human, \file{mgi} for mouse or \file{sgd} for yeast.
#'
#' @return A tibble with columns \code{gene_symbol}, \code{uniprot_id} and \code{term_id}.
fetch_go_genes_go <- function(species) {
  # Binding variables from non-standard evaluation locally
  gene_synonym <- db_object_synonym <- gene_symbol <- symbol <- NULL
  uniprot_id <- db_id <- term_id <- go_term <- NULL


  gaf_file <- stringr::str_glue("http://current.geneontology.org/annotations/{species}.gaf.gz")
  if (!url_exists(gaf_file))
    stop(stringr::str_glue("Cannot download {gaf_file}. Check if {species} is correct species name. Or, perhaps, the server is down. Or your internet is not working."))

  readr::read_tsv(gaf_file, comment = "!", quote = "", col_names = GAF_COLUMNS, show_col_types = FALSE) |>
    dplyr::mutate(gene_synonym = stringr::str_remove(db_object_synonym, "\\|.*$")) |>
    dplyr::select(gene_symbol = symbol, gene_synonym, uniprot_id = db_id, term_id = go_term) |>
    dplyr::distinct()
}


#' Get functional term data from gene ontology
#'
#' Download term information (GO term ID and name) and gene-term mapping
#' (gene symbol and GO term ID) from gene ontology.
#'
#' @param species Species designation. Base file name for species file under
#'   \url{http://current.geneontology.org/annotations}. Examples are
#'   \file{goa_human} for human, \file{mgi} for mouse or \file{sgd} for yeast.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#'
#' @examples
#' \dontrun{
#' go_data <- fetch_go_from_go("sgd")
#' }
fetch_go_from_go <- function(species) {
  terms <- fetch_go_terms()
  mapping <- fetch_go_genes_go(species)

  list(
    terms = terms,
    mapping = mapping
  )
}



#' Download GO term gene mapping from Ensembl
#'
#' @param mart Object class \code{Mart} representing connection to BioMart
#'   database, created with, e.g., \code{useEnsembl}.
#'
#' @return A tibble with columns \code{ensembl_gene_id}, \code{gene_symbol} and
#'   \code{term_id}.
fetch_go_genes_bm <- function(mart) {
  # Binding variables from non-standard evaluation locally
  gene_symbol <- external_gene_name <- term_id <- go_id <- NULL

  biomaRt::getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "go_id"),
    mart = mart
  ) |>
    dplyr::rename(
      gene_symbol = external_gene_name,
      term_id = go_id
    ) |>
    dplyr::filter(term_id != "") |>
    tibble::as_tibble()
}



#' Get functional term data from Ensembl
#'
#' Download term information (GO term ID and name) and gene-term mapping
#' (gene ID, symbol and GO term ID) from Ensembl.
#'
#' @param mart Object class \code{Mart} representing connection to BioMart
#'   database, created with, e.g., \code{useEnsembl}.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @import assertthat
#'
#' @examples
#' \dontrun{
#' mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#' go_terms <- fetch_go_from_bm(mart)
#' }
fetch_go_from_bm <- function(mart) {
  assert_that(class(mart) == "Mart")

  terms <- fetch_go_terms()
  mapping <- fetch_go_genes_bm(mart)

  list(
    terms = terms,
    mapping = mapping
  )
}

