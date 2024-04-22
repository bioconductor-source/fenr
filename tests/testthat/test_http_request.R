timeout_url <- "http://www.google.com:81"


test_that("Correct reaction to a time-out", {
  resp <- http_request(timeout_url, "", timeout = 1)
  expect_is(resp, "list")
  expect_equal(resp$is_error, TRUE)
  expect_equal(resp$status, "unknown")
})
