context("max-min-normalization")

test_that("test normalize sanities", {
  expect_error(normalizeMaxMin(letters[1:10]))
})

test_that("length is constant", {
  x <- rnorm(4L)
  expect_equal(length(x), length(normalizeMaxMin(x)))
})

