library(testthat)

# Set 100 features
N <- 100
features_all <- sprintf("gene_%03d", 1:N)

# Three terms with 50, 10 and 3 features
term_ids <- c("term_1", "term_2", "term_3")
term_names <- c("name_1", "name_2", "name_3")
term_sizes <- c(50, 10, 3)


terms <- tibble::tibble(
  term_id = term_ids,
  term_name = term_names
)

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




test_that("Expected normal output", {
  td <- prepare_for_enrichment(terms, mapping, feature_name = "feature_id")

  # Order is not mandatory, so need to sort before comparison
  expect_equal(sort(term_data$term2name), sort(td$term2name))

  p1 <- purrr::map2(td$term2feature, term_data$term2feature, function(f1, f2) {
    expect_equal(sort(f1), sort(f2))
  })

  p2 <- purrr::map2(td$feature2term, term_data$feature2term, function(f1, f2) {
    expect_equal(sort(f1), sort(f2))
  })
})


test_that("Wrong columns in terms", {
  trms <- terms |>
    dplyr::rename(gobble = term_id)

  expect_error(
    prepare_for_enrichment(trms, mapping, feature_name = "feature_id")
  )
})


test_that("Missing column in mapping", {
  mp <- mapping |>
    dplyr::rename(mobble = term_id)

  expect_error(
    prepare_for_enrichment(terms, mp, feature_name = "feature_id")
  )
})


test_that("Wrong feature name", {
  expect_error(
    prepare_for_enrichment(terms, mapping, feature_name = "trouble")
  )
})


test_that("No overlap between all_features and mapping", {
  expect_error(
    prepare_for_enrichment(terms, mapping, all_features = c("a", "b", "c"),  feature_name = "feature_id")
  )
})

