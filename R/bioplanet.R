#' Get functional term data from BioPlanet
#'
#' Download term information (term ID and name) and gene-pathway mapping
#' (NCBI gene ID, gene symbol and pathway ID) from BioPlanet.
#'
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @examples
#' \dontrun{
#' bioplanet_data <- fetch_bp()
#' }
fetch_bp <- function(use_cache = TRUE) {
  # Binding variables from non-standard evaluation locally
  PATHWAY_ID <- PATHWAY_NAME <- GENE_ID <- GENE_SYMBOL <- NULL

  bp_file <- "https://tripod.nih.gov/bioplanet/download/pathway.csv"

  assert_url_path(bp_file)
  lpath <- cached_url_path("bioplanet", bp_file, use_cache)
  paths <- readr::read_csv(lpath, show_col_types = FALSE)

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
