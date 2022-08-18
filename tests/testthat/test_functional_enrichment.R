library(testthat)

# Set 100 features
N <- 100
features_all <- sprintf("gene_%03d", 1:N)

# Three terms with 50, 10 and 3 features
term_ids <- c("term_1", "term_2", "term_3")
term_names <- c("name_1", "name_2", "name_3")
term_sizes <- c(50, 10, 3)

# Prepare data for enrichment
term2name <- term_names |>
  purrr::set_names(term_ids)

# random selection of features for terms
set.seed(666)
mapping <- purrr::map2_dfr(term_ids, term_sizes, function(tid, n) {
  tibble::tibble(
    term_id = tid,
    feature_id = sample(features_all, n)
  )
})

# Feature to terms conversion list
feature2term <- mapping |>
  dplyr::group_by(feature_id) |>
  dplyr::summarise(terms = list(term_id)) |>
  tibble::deframe()

# Term to feature conversion list
term2feature <- mapping |>
  dplyr::group_by(term_id) |>
  dplyr::summarise(features = list(feature_id)) |>
  tibble::deframe()

# final structure required by functional_enrichment
term_data <- list(
  term2name = term2name,
  feature2term = feature2term,
  term2feature = term2feature
) |>
  structure(class = "fterms")

term_stats <- function(tid, N_sel, n_with_sel) {
  N_with <- length(term2feature[[tid]])
  N_without <- N - N_with
  n_without_sel <- N_sel - n_with_sel
  n_expect <- N_with * N_sel / N
  enrichment <- n_with_sel / n_expect

  p <- 1 - stats::phyper(n_with_sel - 1, N_with, N_without, N_sel)

  sel_with <- term2feature[[tid]][1:n_with_sel]
  if (n_without_sel > 0) {
    sel_without <- setdiff(features_all, term2feature[[tid]])[1:n_without_sel]
  } else {
    sel_without = character(0)
  }
  sel <- c(sel_with, sel_without)

  list(
    N_with = N_with,
    n_with_sel = n_with_sel,
    n_expect = n_expect,
    enrichment = enrichment |> signif(3),
    p = p |> signif(3),
    sel = sel
  )
}



test_that("Enrichment of 1", {
  tid <- term_ids[1]
  N_sel <- 20
  n_with_sel <- 10

  res <- term_stats(tid, N_sel, n_with_sel)
  enr <- functional_enrichment(features_all, res$sel, term_data, fdr_limit = 1) |>
    dplyr::filter(term_id == tid)

  expect_equal(
    c(res$N_with, res$n_with_sel, res$n_expect, res$enrichment, res$p),
    c(enr$N_with, enr$n_with_sel, enr$n_expect, enr$enrichment, enr$p_value)
  )
})


test_that("All selected features in term", {
  tid <- term_ids[2]
  N_sel <- term_sizes[2]
  n_with_sel <- term_sizes[2]

  res <- term_stats(tid, N_sel, n_with_sel)
  enr <- functional_enrichment(features_all, res$sel, term_data, fdr_limit = 1) |>
    dplyr::filter(term_id == tid)

  expect_equal(
    c(res$N_with, res$n_with_sel, res$n_expect, res$enrichment, res$p),
    c(enr$N_with, enr$n_with_sel, enr$n_expect, enr$enrichment, enr$p_value)
  )
})


test_that("No selected features in term", {
  tid <- "term_3"
  N_sel <- 20
  n_with_sel <- 0

  res <- term_stats(tid, N_sel, n_with_sel)
  enr <- functional_enrichment(features_all, res$sel, term_data, fdr_limit = 1) |>
    dplyr::filter(term_id == tid)

  expect_equal(nrow(enr), 0)
})


test_that("Selection is all features", {
  enr <- functional_enrichment(features_all, features_all, term_data, fdr_limit = 1)

  n <- nrow(enr)
  expect_equal(enr$enrichment, rep(1, n))
  expect_equal(enr$p_value, rep(1, n))
})


test_that("No match between selection and all features", {
  expect_null(
    functional_enrichment(features_all, c("a", "b", "c"), term_data)
  )
})


test_that("Incorrect term_data class", {
  td <- term_data
  class(td) <- "wobble"

  expect_error(
    functional_enrichment(features_all, features_all[1:10], td)
  )
})
