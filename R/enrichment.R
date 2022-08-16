#' Prepare term data for enrichment analysis
#'
#' @details Takes simple data frames with functional term information and gene
#'   mapping and converts it into an object required by `functional_enrichment`
#'   for fast analysis.
#'
#' @param term_info Information about term names/descriptions. A tibble with
#'   `term_id` and `term_name`
#' @param mapping Information about term-feature mapping. A tibble with
#'   `term_id` and a feature id, as identified with `feature_name` argument. For
#'   example, if this tibble contains `gene_name` and `term_id`, then
#'   `feature_name = "gene_name"`.
#' @param all_features A vector with all feature ids (background for enrichment).
#' @param feature_name Which column to use from mapping table, e.g. "gene_name" or "gene_id".
#'
#' @return A list required by `functional_enrichment`.
#' @export
prepare_for_enrichment <- function(term_info, mapping, all_features, feature_name = "gene_id") {
  # Check for column name
  if (!(feature_name %in% colnames(mapping))) {
    stop(paste(feature_name, "column not found in mapping table. Check feature_name argument."))
  }

  # Check for missing term descriptions
  mis_term <- setdiff(mapping$term_id, term_info$term_id)
  if (length(mis_term) > 0) {
    dummy <- tibble::tibble(
      term_id = mis_term,
      term_name = rep(NA_character_, length(mis_term))
    )
    term_info <- dplyr::bind_rows(term_info, dummy)
  }

  # List to select term name
  term2name <- term_info$term_name |>
    purrr::set_names(term_info$term_id)

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
  )
}


#' Fast functional enrichment
#'
#' Fast functional enrichment based on hypergeometric distribution. Can be used in interactive applications.
#'
#' @details
#'
#' Functional enrichment in a selection (e.g. significantly DE features) of
#' features, using hypergeometric probability. A feature can be a gene, protein,
#' etc. `term_data` is an object with functional term information and
#' feature-term assignment. It is a list of: `term2info` - a named vector term
#' id => term name; `term2feature` - a list term_id => vector of feature_ids;
#' `feature2term` - a list feature id => vector of term ids. It can be created
#' by `prepare_for_enrichment` function.
#'
#' @param feat_all A character vector with all feature identifiers. This is the
#'   background for enrichment.
#' @param feat_sel A character vector with feature identifiers in the selection.
#' @param term_data Functional term data, as explained in details. It can be
#'   created using \code{prepare_for_enrichment}.
#' @param feat2name An optional named list to convert feature id into feature
#'   name.
#' @param min_count Minimal count of features with term in the selection to be
#'   used.
#' @param fdr_limit Only terms with p_adjust below this limit are returned.
#'
#' @return A tibble with enrichment results. For each term the following
#'   quantities are reported: N_with - number of features with this term in the
#'   among all features, n_with_sel - number of features with this term in the
#'   selection, n_expect - expected number of features with this term in the
#'   selection, under the null hypothesis that terms are assigned to features
#'   randomly, enrichment - ratio of n_with_sel / n_expect; odds_ratio - odds
#'   ratio for enrichment, p_value - p-value from a single hypergeometric test;
#'   p_adjust - p-value adjusted for multiple tests using Benjamini-Hochberg
#'   approach.
#'
#' @examples
#' bp <- fetch_bp()
#' bp_terms <- prepare_for_enrichment(bp$terms, bp$mapping, exmpl_all, feature_name = "gene_name")
#' enr <- functional_enrichment(exmpl_all, exmpl_sel, bp_terms)
#'
#' @export
functional_enrichment <- function(feat_all, feat_sel, term_data, feat2name = NULL,
                                  min_count = 2, fdr_limit = 0.05) {

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
        across(c(N_with, n_with_sel), as.integer),
        across(c(n_expect, enrichment, odds_ratio, p_value), as.numeric),
        p_adjust = p.adjust(p_value, method = "BH"),
        across(c(enrichment, odds_ratio, p_value, p_adjust), ~signif(.x, 3)),
        n_expect = round(n_expect, 2)
      ) |>
      dplyr::arrange(desc(odds_ratio)) |>
      dplyr::filter(p_adjust <= fdr_limit)
  } else {
    NULL
  }
}


