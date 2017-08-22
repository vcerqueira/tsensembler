context("ade basics")

specs <- model_specs(learner = c("bm_svr", "bm_glm", "bm_mars"),
                     learner_pars = NULL)

data("water_consumption")
dataset <- embed_timeseries(water_consumption, 5)
train <- dataset[1:1000,]
validation <- dataset[1001:1200,]
test <- dataset[1201:1500,]

capture.output(model <- ADE(target ~ ., train, specs))

test_that("class ADE", {
  expect_s4_class(model, "ADE")
})

test_that("class base ens", {
  expect_s4_class(model@base_ensemble, "base_ensemble")
})

test_that("class metamodel", {
  expect_true(is.list(model@meta_model))
})

test_that("class metamodel", {
  expect_true(model@omega <= 1 & model@omega > 0)
})

test_that("failing omega", {
  expect_error(ADE(target ~ ., train, specs, omega = 1.1))
})

test_that("failing lambda 1", {
  expect_error(ADE(target ~ ., train, specs, lambda = -1))
})

test_that("failing lambda 2", {
  expect_error(ADE(target ~ ., train, specs, lambda = nrow(train)+1))
})
