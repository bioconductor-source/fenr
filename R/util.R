HTTP_OK <- 200
HTTP_FOUND <- 302


#' Check if URL exists
#'
#' Stops with error message if URL is not accessible.
#'
#' @param url URL to test
#'
#' @return Nothing
assert_url <- function(url) {
  if (!RCurl::url.exists(url))
    stop(stringr::str_glue("URL {url} cannot be found."))
}

#' Check if HTTP file exists
#'
#' Stops with error message if HTTP file is not accessible.
#'
#' @param url_file HTTP address of the file
#'
#' @return Nothing
assert_http_file <- function(url_file) {
  hd <- httr::HEAD(url_file)
  status <- hd$all_headers[[1]]$status
  if (!(status %in% c(HTTP_OK, HTTP_FOUND)))
    stop(stringr::str_glue("HTTP file {url_file} cannot be found. Status = {status}."))
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
