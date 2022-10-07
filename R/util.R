HTTP_OK <- 200
HTTP_FOUND <- 302


#' Check if a tibble contains all required columns
#'
#' @param tb Tibble or data frame to examine
#' @param cols A string vector with required column names
#'
#' @return TRUE if assertion passed
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
#' @return TRUE if assertion passed
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
#' @param fetch_fun A function to retrieve available species, must return a
#'   tibble with a column \code{designation}.
#'
#' @return A tibble with valid species - a response from \code{fetch_fun}
assert_species <- function(species, fetch_fun) {
  assert_that(is.string(species))
  valid_species <- fetch_fun()
  fun_name <- deparse(substitute(fetch_fun))
  assert_that(species %in% valid_species$designation,
    msg = stringr::str_glue("Invalid species {species}. Use {fun_name}() to find all available species.")
  )
  valid_species
}


#' Return a field for a given species
#'
#' Used to match species designation with a species ID.
#'
#' @param species Species designation
#' @param fetch_fun  A function to retrieve available species, must return a
#'   tibble with a column \code{designation}.
#' @param col_name Column name in the tibble returned by \code{fetch_fun} to
#'   extract, e.g. \code{tax_id}
#'
#' @return A value extracted from column \code{col_name} at row where
#'   \code{designation} = \code{species}.
match_species <- function(species, fetch_fun, col_name) {
  designation <- NULL

  sp <- assert_species(species, fetch_fun)
  sp |>
    dplyr::filter(designation == species) |>
    dplyr::pull(get(col_name))
}
