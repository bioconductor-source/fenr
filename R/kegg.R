#' Find all species available from KEGG
#'
#' @return A tibble, in which column \code{organism} contains species
#'   designations used in function \code{fetch_kegg}.
#' @export
#' @examples
#' spe <- fetch_kegg_species()
fetch_kegg_species <- function() {
  KEGGREST::keggList("organism") |>
    tibble::as_tibble()
}


#' Get functional term data from KEGG
#'
#' Download information (pathway ID and name) and gene-pathway mapping (entrez
#' gene ID, gene symbol and pathway ID) from KEGG. Gene symbols are extracted
#' from gene descriptions. For some species (e.g. yeast), gene symbols are
#' returned instead of entrez IDs and not in gene description. This function is
#' based on BioConductor package \pkg{KEGGREST}.
#'
#' @param species KEGG species code, for example "hsa" for human. The
#'   full list of available KEGG species can be found by using \code{fetch_kegg_species}.
#' @param batch_size Number of pathways sent to KEGG database in one query. The
#'   maximum allowed is 10.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @import assertthat
#'
#' @examples
#' kegg_data <- fetch_kegg("sce")
fetch_kegg <- function(species, batch_size = 10) {
  # Binding variables from non-standard evaluation locally
  path_id <- NULL

  assert_that(is.string(species))
  assert_that(is.count(batch_size))
  assert_that(batch_size <= 10, msg = "batch_size needs to be between 1 and 10")

  lst <- tryCatch(
    KEGGREST::keggList("pathway", species),
    error = function(err)
      stop(stringr::str_glue("There is a problem retrieving KEGG pathways for species '{species}'."))
  )

  terms <- tibble::tibble(
    term_id = names(lst) |> stringr::str_remove("path:"),
    term_name = lst
  )
  pids <- terms$term_id
  batches <- split(pids, ceiling(seq_along(pids) / batch_size))

  pb <- progress::progress_bar$new(total = length(batches))
  mapping <- purrr::map_dfr(batches, function(batch) {
    pws <- tryCatch(
      KEGGREST::keggGet(batch),
      error = function(err)
        stop(stringr::str_glue("There is a problem retrieving KEGG batch':\n{err}"))
    )
    pb$tick()
    purrr::map_dfr(pws, function(pw) {
      if (!is.null(pw$GENE)) {
        # KEGG list of genes is a vector with alternate ID and gene description
        gns <-  pw$GENE
        ids <- gns[seq(1, length(gns) - 1, 2)]
        genes <- gns[seq(2, length(gns), 2)] |>
          stringr::str_remove(";.*$")  # attempt to extract gene name
        tibble::tibble(gene_id = ids, gene_symbol = genes, term_id = unname(pw$ENTRY))
      }
    })
  })

  list(
    terms = terms,
    mapping = mapping
  )
}

