#' Prepare term data for enrichment analysis
#'
#' Prepare term data downloaded with \code{fetch_*} functions for fast
#' enrichment analysis.
#'
#' @details
#'
#' Takes two data frames with functional term information (\code{terms}) and
#' gene mapping (\code{mapping}) and converts them into an object required by
#' \code{functional_enrichment} for fast analysis. Terms and mapping can be
#' created with database access functions in this package, for example
#' \code{fetch_reactome} or \code{fetch_go_from_go}.
#'
#' @param terms Information about term names/descriptions. A tibble with columns
#'   \code{term_id} and \code{term_name}.
#' @param mapping Information about term-feature mapping. A tibble with
#'   \code{term_id} and a feature id, as identified with \code{feature_name}
#'   argument. For example, if this tibble contains \code{gene_symbol} and
#'   \code{term_id}, then you need to set \code{feature_name = "gene_symbol"}.
#' @param all_features A vector with all feature ids used as background for
#'   enrichment. If not specified, all features found in \code{mapping} will be
#'   used, resulting in a larger object size.
#' @param feature_name Which column to use from mapping table, e.g.
#'   \code{gene_symbol} or \code{ensembl_gene_id}.
#'
#' @return An object class \code{fenr_terms} required by
#'   \code{functional_enrichment}.
#' @export
#'
#' @examples
#' bp <- fetch_bp()
#' bp_terms <- prepare_for_enrichment(bp$terms, bp$mapping, exmpl_all, feature_name = "gene_symbol")
prepare_for_enrichment <- function(terms, mapping, all_features = NULL, feature_name = "gene_id") {
  # Binding variables from non-standard evaluation locally
  feature_id <- term_id <- NULL

  # Check terms
  if (!all(c("term_id", "term_name") %in% colnames(terms)))
    stop("Column names in 'terms' should be 'term_id' and 'term_name'.")

  # Check mapping
  if (!("term_id" %in% colnames(mapping)))
    stop("'mapping' should contain a column named 'term_id'.")

  # Check for feature name
  if (!(feature_name %in% colnames(mapping)))
    stop(feature_name, "column not found in mapping table. Check feature_name argument.")

  # Replace empty all_features with everything from mapping
  map_features <- mapping[[feature_name]] |>
    unique()
  if (is.null(all_features)) {
    all_features <- map_features
  } else {
    # Check if mapping is contained in all features
    if (length(intersect(all_features, map_features)) == 0)
      stop("No overlap between 'all_features' and features found in 'mapping'. Did you provide correct 'all_features'?")
  }

  # Check for missing term descriptions
  mis_term <- setdiff(mapping$term_id, terms$term_id)
  if (length(mis_term) > 0) {
    dummy <- tibble::tibble(
      term_id = mis_term,
      term_name = rep(NA_character_, length(mis_term))
    )
    terms <- dplyr::bind_rows(terms, dummy)
  }

  # List to select term name
  term2name <- terms$term_name |>
    purrr::set_names(terms$term_id)

  # feature-term tibble
  feature_term <- mapping |>
    dplyr::rename(feature_id = !!feature_name) |>
    dplyr::filter(feature_id %in% all_features)

  # Feature to terms conversion list
  feature2term <- feature_term |>
    dplyr::group_by(feature_id) |>
    dplyr::summarise(terms = list(term_id)) |>
    tibble::deframe()

  # Term to feature conversion list
  term2feature <- feature_term |>
    dplyr::group_by(term_id) |>
    dplyr::summarise(features = list(feature_id)) |>
    tibble::deframe()

  list(
    term2name = term2name,
    term2feature = term2feature,
    feature2term = feature2term
  ) |>
    structure(class = "fenr_terms")
}


#' Fast functional enrichment
#'
#' Fast functional enrichment based on hypergeometric distribution. Can be used
#' in interactive applications.
#'
#' @details
#'
#' Functional enrichment in a selection (e.g. significantly DE features) of
#' features, using hypergeometric probability (that is Fisher's exact test). A
#' feature can be a gene, protein, etc. \code{term_data} is an object with
#' functional term information and feature-term assignment.
#'
#' @param feat_all A character vector with all feature identifiers. This is the
#'   background for enrichment.
#' @param feat_sel A character vector with feature identifiers in the selection.
#' @param term_data An object class \code{fenr_terms}, created by
#'   \code{prepare_for_enrichment}. It is a list of three elements: \itemize{
#'   \item{\code{term2info} - a named vector term id => term name}
#'   \item{\code{term2feature} - a list term id => vector of feature ids}
#'   \item{\code{feature2term} - a list feature id => vector of term ids}}
#' @param feat2name An optional named list to convert feature id into feature
#'   name.
#' @param min_count Minimal count of features with term in the selection to be
#'   used.
#' @param fdr_limit Only terms with p_adjust below this limit are returned.
#'
#' @return A tibble with enrichment results. For each term the following
#'   quantities are reported: \itemize{ \item{\code{N_with} - number of features
#'   with this term in the among all features} \item{\code{n_with_sel} - number of
#'   features with this term in the selection} \item{\code{n_expect} - expected
#'   number of features with this term in the selection, under the null
#'   hypothesis that terms are assigned to features randomly}
#'   \item{\code{enrichment} - ratio of n_with_sel / n_expect}
#'   \item{\code{odds_ratio} - odds ratio for enrichment; is infinite, when all
#'   features with the given term are in the selection} \item{\code{p_value} -
#'   p-value from a single hypergeometric test} \item{\code{p_adjust} - p-value
#'   adjusted for multiple tests using Benjamini-Hochberg approach} }
#'
#' @examples
#' bp <- fetch_bp()
#' bp_terms <- prepare_for_enrichment(bp$terms, bp$mapping, exmpl_all, feature_name = "gene_symbol")
#' enr <- functional_enrichment(exmpl_all, exmpl_sel, bp_terms)
#'
#' @import assertthat
#' @importFrom methods is
#' @export
functional_enrichment <- function(feat_all, feat_sel, term_data, feat2name = NULL,
                                  min_count = 2, fdr_limit = 0.05) {

  # Binding variables from non-standard evaluation locally
  N_with <- n_with_sel <- n_expect <- enrichment <- odds_ratio <- NULL
  desc <- p_value <- p_adjust <- NULL

  # Check term_data class
  assert_that(is(term_data, "fenr_terms"))

  # If no overlap between selection and all, return NULL
  if(!any(feat_sel %in% feat_all))
    return(NULL)

  # all terms present in the selection
  our_terms <- feat_sel |>
    purrr::map(\(x) term_data$feature2term[[x]]) |>
    unlist() |>
    unique()
  # number of features in selection
  N_sel <- length(feat_sel)
  # total number of features
  N_tot <- length(feat_all)

  res <- purrr::map_dfr(our_terms, function(term_id) {
    # all features with the term
    tfeats <- term_data$term2feature[[term_id]]

    # features from selection with the term
    tfeats_sel <- intersect(tfeats, feat_sel)

    N_with <- length(tfeats)
    N_without <- N_tot - N_with

    # building contingency table
    n_with_sel <- length(tfeats_sel)
    n_without_sel <- N_sel - n_with_sel
    n_with_nsel <- N_with - n_with_sel
    n_without_nsel <- N_tot - (n_with_sel + n_without_sel + n_with_nsel)

    # contingency table is
    #               | In selection  | Not in selection
    #-------------------------------------------------
    #  With term    | n_with_sel    | n_with_nsel
    #  Without term | n_without_sel | n_without_nsel

    if (n_with_sel < min_count) return(NULL)

    # Expected number of features in selection, if random
    n_expect <- N_with * N_sel / N_tot

    # Odds ratio
    odds_ratio <- (n_with_sel / n_without_sel) / (n_with_nsel / n_without_nsel)

    # Hypergeometric function much faster than fisher.test
    p <- 1 - stats::phyper(n_with_sel - 1, N_with, N_without, N_sel)

    if (!is.null(feat2name)) tfeats_sel <- feat2name[tfeats_sel] |> unname()

    term_name <- term_data$term2name[[term_id]]
    # returns NAs if no term found
    if (is.null(term_name)) term_name <- NA_character_

    # constructing a tibble in every iteration is more time expensive than `c`,
    # even with overhead of converting types afterwards. Not elegant, but fast.
    c(
      term_id = term_id,
      term_name = term_name,
      N_with = N_with,
      n_with_sel = n_with_sel,
      n_expect = n_expect,
      enrichment = n_with_sel / n_expect,
      odds_ratio = odds_ratio,
      ids = paste(tfeats_sel, collapse = ", "),
      p_value = p
    )
  })
  # Drawback - if all selections below minimum, res is tibble 0 x 0, need to catch it
  if (nrow(res > 0)) {
    res |>
      dplyr::mutate(
        dplyr::across(c(N_with, n_with_sel), as.integer),
        dplyr::across(c(n_expect, enrichment, odds_ratio, p_value), as.numeric),
        p_adjust = stats::p.adjust(p_value, method = "BH"),
        dplyr::across(c(enrichment, odds_ratio, p_value, p_adjust), ~signif(.x, 3)),
        n_expect = round(n_expect, 2)
      ) |>
      dplyr::arrange(desc(odds_ratio)) |>
      dplyr::filter(p_adjust <= fdr_limit)
  } else {
    NULL
  }
}


