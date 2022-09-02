HTTP_OK <- 200
HTTP_FOUND <- 302

#' Check if URL with a path is responding OK
#'
#' Stops with error message if the path is not accessible or server not
#' responding.
#'
#' @param url_path Full URL with a path, e.g. `https://reactome.org/download/current/ReactomePathways.txt`.
#'
#' @return Nothing
assert_url_path <- function(url_path) {
  assert_that(is.string(url_path))
  hd <- tryCatch(
    httr::HEAD(url_path),
    error = function(e) {
      stop(e)
    }
  )
  status <- hd$all_headers[[1]]$status
  if (!(status %in% c(HTTP_OK, HTTP_FOUND)))
    stop(stringr::str_glue("HTTP path {url_path} cannot be found. Status = {status}."))
}



#' Check if species are valid
#'
#' @param species A string, species designation for a given database
#' @param fetch_fun A function to retrieve available species, must return a
#'   tibble with a column \code{designation}.
#'
#' @return TRUE is correct
assert_species <- function(species, fetch_fun) {
  assert_that(is.string(species))
  valid_species <- fetch_fun()
  fun_name <- deparse(substitute(fetch_fun))
  assert_that(species %in% valid_species$designation,
    msg = stringr::str_glue("Invalid species {species}. Use {fun_name}() to find all available species.")
  )
}
