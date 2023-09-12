#' Retrieve Cache Location for `fenr`
#'
#' This function returns the cache directory location for the `fenr` package.
#'
#' @return A character string representing the path to the cache directory for `fenr`.
#' @importFrom tools R_user_dir
cache_location <- function() {
  tools::R_user_dir("fenr", which = "cache")
}


#' URL Cache Path Function
#'
#' This function queries a BiocFileCache object for a given resource name (`rname`).
#' If the resource isn't found, it's added to the cache. If there's one entry found,
#' it checks for updates and retrieves the local path. If multiple entries are found,
#' it cleans up by removing them and then adds the resource anew.
#'
#' @param bfc A BiocFileCache object.
#' @param rname A character string representing the resource name to be queried or added.
#' @param fpath A character string representing the path to the file.
#'
#' @return A character string representing the local path of the cached resource.
url_cache_path <- function(bfc, rname, fpath) {
  hits <- BiocFileCache::bfcquery(bfc, rname, field = "rname")

  # Not found, add to cache, get local path
  if(nrow(hits) == 0) {
    lpath <- BiocFileCache::bfcadd(bfc, rname, fpath)
  }

  # Found one entry, check for updates, get local path
  else if(nrow(hits) == 1) {
    if(BiocFileCache::bfcneedsupdate(bfc, hits$rid))
      BiocFileCache::bfcupdate(bfc, hits$rid, fpath = fpath, ask = FALSE)
    lpath <- BiocFileCache::bfcpath(bfc, hits$rid)
  }

  # multiple entries found, clean up
  else {
    BiocFileCache::bfcremove(bfc, hits$rid)
    lpath <- BiocFileCache::bfcadd(bfc, rname, fpath)
  }

  return(lpath)
}


#' Cached URL Path Function
#'
#' This function allows you to obtain a cached local path for a resource identified
#' by its name and file path. If caching is enabled (controlled by the `use_cache`
#' argument), it uses a BiocFileCache object to store and retrieve the resource.
#' If caching is disabled, the original file path is returned.
#'
#' @param rname A character string representing the resource name.
#' @param fpath A character string representing the path to the file.
#' @param use_cache A logical value indicating whether to use caching. If TRUE,
#'   the function caches the resource; if FALSE, it uses the original path.
#'
#' @return A character string representing either the cached local path (if caching
#'   is enabled) or the original file path (if caching is disabled).
cached_url_path <- function(rname, fpath, use_cache) {
  if(use_cache) {
    cache <- fenr_cache_location()
    bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
    lpath <- url_cache_path(bfc, rname, fpath)
  } else {
    lpath <- fpath
  }

  return(lpath)
}
