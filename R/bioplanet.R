#' URL of Bioplanet top dir
#'
#' @return A string with URL.
#' @noRd
get_bioplanet_download <- function() {
  getOption("BIOPLANET_PATHWAY_DOWNLOAD", "https://tripod.nih.gov/bioplanet/download")
}

#' URL of Bioplanet pathway file
#'
#' @return A string with URL.
#' @noRd
get_bioplanet_pathway_file <- function() {
  getOption("BIOPLANET_PATHWAY_FILE", "https://tripod.nih.gov/bioplanet/download/pathway.csv")
}


#' Get functional term data from BioPlanet
#'
#' Download term information (term ID and name) and gene-pathway mapping
#' (NCBI gene ID, gene symbol and pathway ID) from BioPlanet.
#'
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @examples
#' bioplanet_data <- fetch_bp(on_error = "warn")
fetch_bp <- function(use_cache = TRUE, on_error = c("stop", "warn", "ignore")) {
  on_error <- match.arg(on_error)

  # Binding variables from non-standard evaluation locally
  PATHWAY_ID <- PATHWAY_NAME <- GENE_ID <- GENE_SYMBOL <- NULL

  if(!assert_url_path(get_bioplanet_download(), on_error))
    return(NULL)

  pathway_file <- get_bioplanet_pathway_file()

  lpath <- cached_url_path("bioplanet_pathway", pathway_file, use_cache)
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
