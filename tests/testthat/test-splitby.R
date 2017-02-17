context("splitby")

test_that("check splitby", {
  .expr <- "welele.welele"
  expect_match(.splitBy.(.expr)[1], "welele")
})

