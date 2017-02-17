context("fifo")

test_that("mode in value equal mode x", {
  x <- rnorm(4L)
  y <- 5.
  expect_equal(mode(x), mode(y))
})

test_that("fifo length is constant", {
  x <- rnorm(4L)
  y <- 5.
  expect_equal(length(x), length(FIFO(x, y)))
})
