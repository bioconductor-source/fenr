#' Get functional term data from BioPlanet
#'
#' Download term information (term ID and name) and gene-pathway mapping
#' (NCBI gene ID, gene symbol and pathway ID) from BioPlanet.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @examples
#' bioplanet_data <- fetch_bp()
fetch_bp <- function() {
  # Binding variables from non-standard evaluation locally
  PATHWAY_ID <- PATHWAY_NAME <- GENE_ID <- GENE_SYMBOL <- NULL

  bp_file <- "https://tripod.nih.gov/bioplanet/download/pathway.csv"

  # ---------------------------
  # WARNING: A NASTY HACK
  #
  # SSL certificate for tripod.nih.gov expired, so normal read_csv does not work.
  # The workaround is to use httr::get in an unsecure SSL disabled environment.
  # res <- httr::with_config(
  #  config = httr::config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE),
  #  httr::GET(bp_file)
  # )
  # paths <- httr::content(res, show_col_types = FALSE, encoding = "UTF-8")

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
