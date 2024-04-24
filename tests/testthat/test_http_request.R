timeout_url <- "http://www.google.com:81"


test_that("Correct reaction to a time-out in http_resuest()", {
  resp <- http_request(timeout_url, "", timeout = 1)
  expect_is(resp, "list")
  expect_equal(resp$is_error, TRUE)
  expect_equal(resp$status, "unknown")
})

test_that("Correct reaction to a time-out in assert_url_path()", {
  expect_warning({assert_url_path(timeout_url, on_error = "warn", timeout = 1)})
  resp <- assert_url_path(timeout_url, on_error = "ignore", timeout = 1)
  expect_equal(resp, FALSE)
})
