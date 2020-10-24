## tsensembler

An S4 package for dynamic combination of forecasting models

### Dynamic Ensembles for Time Series Forecasting

The package provided methods for dynamically combining forecasting models for time series forecasting predictive tasks. It leverages machine learning models from other packages to automatically combine expert advice using metalearning and other state-of-the-art forecasting combination approaches. 

## Installing

Install the package using your R console:

`install.packages('tsensembler')` 

or

`devtools::install_github('vcerqueira/tsensembler')`


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

nall_kernels <- c("rbfdot","polydot","vanilladot","laplacedot")

pars_predictors <- list(
  bm_gaussianprocess = list(kernel = nall_kernels, tol = c(.001)),
  bm_svr = list(kernel = nall_kernels, C=c(1, 5), epsilon=c(.1,0.01)),
  bm_ppr = list(nterms = c(2,4),
                sm.method = c("supsmu","gcvspline")),
  bm_mars = list(degree = c(1, 3), nk = c(10,20),
                 thresh=c(0.001),
                 pmethod=c("forward")),
  bm_glm = list(alpha = c(0,.25,.5,.75,1),
                family = c("gaussian")),
  bm_randomforest = list(num.trees = c(250,500),
                         mtry = c(5,10)),
  bm_pls_pcr = list(method = c("simpls","kernelpls","svdpc")),
  bm_cubist  = list(committees= c(1,5, 15)),
  bm_xgb = list()
)


base_predictors <- names(pars_predictors)

specs <- model_specs(base_predictors,pars_predictors)

# building the ensemble
model <- quickADE(target ~., train, specs)

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

Any bug report or suggestion please contact me at cerqueira.vitormanuel@gmail.com
