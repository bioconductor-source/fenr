HTTP_OK <- 200
HTTP_FOUND <- 302

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

GAF_TYPES <- rep("c", length(GAF_COLUMNS)) |>
  stringr::str_c(collapse = "")


#' Check if a tibble contains all required columns
#'
#' @param tb Tibble or data frame to examine
#' @param cols A string vector with required column names
#'
#' @return TRUE if assertion passed
#' @noRd
assert_columns <- function(tb, cols) {
  if(!all(cols %in% colnames(tb)))
    stop(paste0("The data frame needs to contain the following columns\n", paste(cols, collapse = ", ")))
  return(TRUE)
}

#' Check if URL with a path is responding OK
#'
#' Stops with error message if the path is not accessible or server not
#' responding.
#'
#' @param url_path Full URL with a path, e.g. `https://reactome.org/download/current/ReactomePathways.txt`.
#' @param stop_if_error Logical, if TRUE stops with error message, otherwise returns FALSE
#'
#' @importFrom assertthat assert_that is.string
#' @return TRUE if assertion passed
#' @noRd
assert_url_path <- function(url_path, stop_if_error = TRUE) {
  assert_that(is.string(url_path))
  hd <- tryCatch(
    httr::HEAD(url_path),
    error = function(e) {
      stop(e)
    }
  )
  status <- hd$all_headers[[1]]$status
  if (!(status %in% c(HTTP_OK, HTTP_FOUND))) {
    if(stop_if_error) {
      stop(stringr::str_glue("HTTP path {url_path} cannot be found. Status = {status}."))
    } else {
      return(FALSE)
    }
  }
  return(TRUE)
}



#' Check if species are valid
#'
#' @param species A string, species designation for a given database
#' @param fetch_fun A string, name of the function to retrieve available
#'   species, must return a tibble with a column \code{designation}.
#'
#' @importFrom assertthat assert_that is.string
#' @return A tibble with valid species - a response from \code{fetch_fun}
#' @noRd
assert_species <- function(species, fetch_fun) {
  assert_that(is.string(species))
  assert_that(is.string(fetch_fun))
  f <- match.fun(fetch_fun)
  valid_species <- f()
  assert_that(species %in% valid_species$designation,
    msg = stringr::str_glue("Invalid species {species}. Use {fetch_fun}() to find all available species.")
  )
  valid_species
}


#' Return a field for a given species
#'
#' Used to match species designation with a species ID.
#'
#' @param species Species designation
#' @param fetch_fun  A string with name of a function to retrieve available
#'   species, must return a tibble with a column \code{designation}.
#' @param col_name Column name in the tibble returned by \code{fetch_fun} to
#'   extract, e.g. \code{tax_id}
#'
#' @return A value extracted from column \code{col_name} at row where
#'   \code{designation} = \code{species}.
#' @noRd
match_species <- function(species, fetch_fun, col_name) {
  designation <- NULL

  sp <- assert_species(species, fetch_fun)
  sp |>
    dplyr::filter(designation == species) |>
    dplyr::pull(get(col_name))
}


#' Get features annotated with a given term.
#'
#' @param term_data An object class \code{fenr_terms}, created by
#'   \code{prepare_for_enrichment}.
#' @param term_id A string with a functional term ID.
#'
#' @return A character vector containing feature IDs annotated with the term ID.
#'
#' @importFrom methods is
#' @importFrom assertthat assert_that is.string
#' @export
#' @examples
#' go_data <- fetch_go(species = "sgd")
#' go_terms <- prepare_for_enrichment(go_data$terms, go_data$mapping, feature = "gene_symbol")
#' feats <- get_term_features(go_terms, "GO:0000001")
get_term_features <- function(term_data, term_id) {
  # Check term_data class
  assert_that(is(term_data, "fenr_terms"))
  assert_that(is.string(term_id))

  term_data$term2feature[[term_id]]
}


#' Get terms annotating a given feature.
#'
#' @param term_data An object class \code{fenr_terms}, created by
#'   \code{prepare_for_enrichment}.
#' @param feature_id A string with a feature ID
#'
#' @return A character vector containing functional term IDs annotating given feature.
#' @importFrom assertthat assert_that is.string
#' @importFrom methods is
#' @export
#' @examples
#' go_data <- fetch_go(species = "sgd")
#' go_terms <- prepare_for_enrichment(go_data$terms, go_data$mapping, feature = "gene_symbol")
#' trms <- get_feature_terms(go_terms, "GEM1")
get_feature_terms <- function(term_data, feature_id) {
  # Check term_data class
  assert_that(is(term_data, "fenr_terms"))
  assert_that(is.string(feature_id))

  term_data$feature2term[[feature_id]]
}



#' Test correctness of data structure returned by fetch_*
#'
#' @param re The data structure returned by fetch_*
#' @noRd
test_fetched_structure <- function(re) {
  expect_is(re, "list")
  expect_length(re, 2)
  expect_named(re)
  expect_equal(names(re), c("terms", "mapping"))
}


#' Test if returned terms tibble is consistent with expected selection
#'
#' @param returned A tibble, the terms data returned by  fetch_*.
#' @param expected A tibble, a selection of expected terms data.
test_terms <- function(returned, expected) {
  expect_is(returned, "tbl")
  expect_contains(names(returned), c("term_id", "term_name"))

  returned <- returned |> dplyr::select(term_id, term_name)

  merged <- expected |>
    dplyr::left_join(returned, by = c("term_id", "term_name")) |>
    tidyr::drop_na()
  expect_equal(nrow(expected), nrow(merged))
}

#' Test if returned mapping is consistent with expected selection
#'
#' @param returned A tibble, the mapping returned by  fetch_*.
#' @param expected A tibble, a selection of expected mapping.
test_mapping <- function(returned, expected, feature_id) {
  expect_is(returned, "tbl")
  expect_contains(names(returned), c("term_id", feature_id))

  returned <- returned |> dplyr::select(term_id, all_of(feature_id))

  merged <- expected |>
    dplyr::left_join(returned, by = c("term_id", feature_id)) |>
    dplyr::distinct() |>
    tidyr::drop_na()
  expect_equal(nrow(expected), nrow(merged))
}


