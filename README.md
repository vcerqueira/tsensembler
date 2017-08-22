## tsensembler

An S4 package for dynamic combination of forecasting models

### Dynamic Ensembles for Time Series Forecasting

The package provided methods for dynamically combining forecasting models for time series forecasting predictive tasks. It leverages machine learning models from other packages to automatically combine expert advice using metalearning and other state-of-the-art forecasting combination approaches. 

## Installing

Install the package using your R console:

`install.packages('tsensembler')`


## Illustrative examples

```
# Using data of water consumption time series attached to the package
data("water_consumption")

embedding time series into a matrix
`dataset <- embed_timeseries(water_consumption, 5)`

# splitting data into train/test
train <- dataset[1:1000,]
test <- dataset[1001:1020, ]

# setting up base model parameters
specs <- model_specs(
  learner = c("bm_ppr","bm_glm","bm_svr","bm_mars"), 
  learner_pars = list(
    bm_glm = list(alpha = c(0, .5, 1)),
    bm_svr = list(kernel = c("rbfdot", "polydot"),
                  C = c(1,3)),
    bm_ppr = list(nterms = 4)
  ))

# building the ensemble
model <- ADE(target ~., train, specs)

# forecast next value and update base and meta models
# every three points;
# in the other points, only the weights are updated
predictions <- numeric(nrow(test))
for (i in seq_along(predictions)) {
  predictions[i] <- predict(model, test[i, ])@y_hat
  if (i %% 3 == 0) {
    model <-
      update_base_models(model,
                         rbind.data.frame(train, test[seq_len(i), ]))

    model <- update_ade_meta(model, rbind.data.frame(train, test[seq_len(i), ]))
  }
  else
    model <- update_weights(model, test[i, ])
}

point_forecast <- forecast(model, h = 5)

# setting up an ensemble of support vector machines
specs2 <-
  model_specs(learner = c("bm_svr"),
              learner_pars = list(
                bm_svr = list(kernel = c("vanilladot", "polydot",
                                         "rbfdot"),
                              C = c(1,3,6))
              ))

model <- DETS(target ~., train, specs2)
preds <- predict(model, test)@y_hat
```

## Contact

Any bug report or suggestion please contact me at vmac@inesctec.pt
