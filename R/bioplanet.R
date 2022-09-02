#' Get functional term data from BioPlanet
#'
#' Download term information (term ID and name) and gene-pathway mapping
#' (NCBI gene ID, gene symbol and pathway ID) from BioPlanet.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#'
#' @examples
#' bioplanet_data <- fetch_bp()
fetch_bp <- function() {
  # Binding variables from non-standard evaluation locally
  PATHWAY_ID <- PATHWAY_NAME <- GENE_ID <- GENE_SYMBOL <- NULL
  term_id <- term_name <- ncbi_id <- gene_symbol <- NULL

  bp_file <- "https://tripod.nih.gov/bioplanet/download/pathway.csv"
  assert_url_path(bp_file)
  paths <- readr::read_csv(bp_file, show_col_types = FALSE)

  terms <- paths |>
    dplyr::select(term_id = PATHWAY_ID, term_name = PATHWAY_NAME) |>
    dplyr::distinct()

  mapping <- paths |>
    dplyr::select(term_id = PATHWAY_ID, ncbi_id = GENE_ID, gene_symbol = GENE_SYMBOL) |>
    dplyr::distinct()

  list(
    terms = terms,
    mapping = mapping
  )
}
