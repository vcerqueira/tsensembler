context("recent series")

specs <- model_specs(learner = c("bm_svr", "bm_glm", "bm_mars"),
                     learner_pars = NULL)

data("water_consumption")
dataset <- embed_timeseries(water_consumption, 5)
train <- dataset[1:1000,]
validation <- dataset[1001:1200,]
test <- dataset[1201:1500,]

capture.output(model <- ADE(target ~ ., train, specs))

RS <- model@recent_series


test_that("colnames recent series", {
  expect_equal(colnames(RS), colnames(train))
})


test_that("dim recent series", {

  expect_equal(nrow(RS), model@lambda)
})
