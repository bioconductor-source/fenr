#' Get functional term data from KEGG
#'
#' Download information (pathway ID and name) and gene-pathway mapping (entrez
#' gene ID, gene symbol and pathway ID) from KEGG. Gene symbols are extracted
#' from gene descriptions. For some species (e.g. yeast), gene symbols are
#' returned instead of entrez IDs and not in gene description. This function is
#' based on BioConductor package \code{KEGGREST}.
#'
#' @param species KEGG species code, for example \code{"hsa"} for human. The
#'   full list of available KEGG species can be found by using
#'   \code{KEGGREST::keggList("organism")}. The column \code{organism} contains
#'   the codes used here.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#'
#' @examples
#' \dontrun{
#' kegg_data <- fetch_kegg("hsa")
#' }
fetch_kegg <- function(species) {
  lst <- KEGGREST::keggList("pathway", species)
  terms <- tibble::tibble(
    term_id = names(lst) |> stringr::str_remove("path:"),
    term_name = lst
  )

  pb <- progress::progress_bar$new(total = nrow(terms))
  mapping <- purrr::map_dfr(terms$term_id, function(path_id) {
    pw <- KEGGREST::keggGet(path_id)
    pb$tick()
    if (!is.null(pw[[1]]$GENE)) {
      # KEGG list of genes is a vector with alternate ID and gene description
      gns <-  pw[[1]]$GENE
      ids <- gns[seq(1, length(gns) - 1, 2)]
      genes <- gns[seq(2, length(gns), 2)] |>
        stringr::str_remove(";.*$")  # attempt to extract gene name
      tibble::tibble(gene_id = ids, gene_symbol = genes, term_id = path_id)
    }
  })

  list(
    terms = terms,
    mapping = mapping
  )
}

