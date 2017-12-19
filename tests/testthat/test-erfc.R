context("erfc function")

test_that("check mode of x", {
  x <- rnorm(4L)
  expect_match(mode(erfc(x)), "numeric")
})

test_that("erfc length is constant", {
  x <- rnorm(4L)
  y <- 5.
  expect_equal(length(x), length(erfc(x)))
})
