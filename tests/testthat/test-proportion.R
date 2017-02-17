context("proportio")

test_that("check mode of x", {
  x <- rnorm(4L)
  expect_equal(sum(proportion(x)), 1)
})

test_that("erfc length is constant", {
  x <- rnorm(4L)
  expect_equal(length(x), length(proportion(x)))
})
