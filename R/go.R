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
  go <- ontologyIndex::get_ontology(obo_file, extract_tags = "minimal")
  tibble::tibble(
    term_id = go$id,
    term_name = go$name
  )
}


#' Download GO term gene mapping
#'
#' @param species Root name for species file under
#'   \code{http://current.geneontology.org/annotations}. Examples are
#'   \code{goa_human} for human, \code{mgi} for mouse or \code{sgd} for yeast.
#'
#' @return A tibble with gene_name, uniprot_id and term_id.
fetch_go_genes_go <- function(species) {
  gaf_file <- stringr::str_glue("http://current.geneontology.org/annotations/{species}.gaf.gz")
  stopifnot(url_exists(gaf_file))

  readr::read_tsv(gaf_file, comment = "!", quote = "", col_names = GAF_COLUMNS, show_col_types = FALSE) |>
    dplyr::mutate(gene_name = stringr::str_remove(db_object_synonym, "\\|.*$")) |>
    dplyr::select(gene_name, uniprot_id = db_id, term_id = go_term) |>
    dplyr::distinct()
}


#' Get functional term data from gene ontology
#'
#' @param species Root name for species file under
#'   \code{http://current.geneontology.org/annotations}. Examples are
#'   \code{"goa_human"} for human, \code{"mgi"} for mouse or \code{"sgd"} for yeast.
#'
#' @return A list with terms and mapping tibbles.
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



fetch_go_genes_bm <- function(mart) {
  gene2go <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "go_id"),
    mart = mart
  ) |>
    rename(
      gene_id = ensembl_gene_id,
      gene_name = external_gene_name,
      term_id = go_id
    ) |>
    filter(term_id != "") |>
    as_tibble()
}


# gene names and ensembl ids
fetch_go_from_bm <- function(mart) {
  terms <- fetch_go_terms()
  mapping <- fetch_go_genes_bm(mart)

  list(
    terms = terms,
    mapping = mapping
  )
}

