context("selectbest")

x <- as.data.frame(matrix(rnorm(20), nrow=5))
xb <- select_best(x)

test_that("ROW SUMS SELECT BEST", {

  expect_equal(unname(rowSums(xb)), rep(1,times=nrow(x)))
})

test_that("output matrix", {

  expect_equal(dim(xb), dim(x))
})
