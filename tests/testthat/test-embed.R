context("embedding")

test_that("check embedding and unembedding", {
  ts1 <- xts::as.xts(rnorm(100L), order.by = Sys.Date() + seq_len(100L))
  ts2 <- xts::as.xts(rnorm(100L), order.by = Sys.Date() + seq_len(100L))
  embedded.ts1 <- embed.timeseries(ts1, 30L)
  embedded.ts2 <- embed.timeseries(ts2, 30L)

  ts.original <- unembed.timeseries(embedded.ts1, embedded.ts2)$train

  expect_equal(length(ts.original), NROW(c(ts1)))
})

