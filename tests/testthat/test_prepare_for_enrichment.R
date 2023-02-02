library(testthat)


# Set 100 features
N <- 100
features_all <- sprintf("gene_%03d", seq_len(N))

# Three terms with 50, 10 and 3 features
term_ids <- c("term_1", "term_2", "term_3")
term_names <- c("name_1", "name_2", "name_3")
term_sizes <- c(50, 10, 3)


terms <- tibble::tibble(
  term_id = term_ids,
  term_name = term_names
)

# random selection of features for terms
set.seed(666)
mapping <- purrr::map2_dfr(term_ids, term_sizes, function(tid, n) {
  tibble::tibble(
    term_id = tid,
    feature_id = sample(features_all, n)
  )
})




test_that("Expected correct output", {
  td <- prepare_for_enrichment(terms, mapping, feature_name = "feature_id")

  # Check term names
  for(i in seq_along(terms$term_id)) {
    r <- terms[i, ]
    expect_equal(sort(r$term_name), sort(td$term2name[[r$term_id]]))
  }


  # Check term-feature hash
  term_ids <- mapping$term_id |> unique()
  chk1 <- term_ids |>
    purrr::map(function(trm) {
      expected <- mapping |>
        dplyr::filter(term_id == trm) |>
        dplyr::pull(feature_id) |>
        sort()
      returned <- td$term2feature[[trm]] |>
        sort()
      expect_equal(expected, returned)
    })

  # Check feature-term hash
  feature_ids <- mapping$feature_id |> unique()
  chk2 <- feature_ids |>
    purrr::map(function(feat) {
      expected <- mapping |>
        dplyr::filter(feature_id == feat) |>
        dplyr::pull(term_id) |>
        sort()
      returned <- td$feature2term[[feat]] |>
        sort()
      expect_equal(expected, returned)
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

