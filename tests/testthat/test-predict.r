context("ade predict workflow")

data("water_consumption")
dataset <- embed_timeseries(water_consumption, 5)

train <- dataset[1:1000,]
test <- dataset[1001:1010, ]

specs <- model_specs(c("bm_ppr","bm_glm","bm_mars"), NULL)
## same with model <- DETS(target ~., train, specs)
capture.output(model <- ADE(target ~., train, specs))

# if consecutive know observations are predicted (e.g. a validation/test set)
# the updating is automatically done internally.
predictions1 <- predict(model, test)@y_hat

# otherwise, the models need to be updated
predictions <- numeric(nrow(test))
# predict new data and update the weights of the model
for (i in seq_along(predictions)) {
  predictions[i] <- predict(model, test[i, ])@y_hat

  model <- update_weights(model, test[i, ])
}


test_that("preds equal", {
  expect_identical(predictions1, predictions)
})

