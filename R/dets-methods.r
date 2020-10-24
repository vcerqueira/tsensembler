setMethod("show",
          signature("DETS"),
          function(object) {
            cat("## Dynamic Ensemble for Time Series Forecasting ##\n\n")

            cat("Base ensemble specs:\n")
            print(object@base_ensemble)
            cat("\n")

            cat("Lambda is set to:", object@lambda, "\n")
            cat("Omega is set to:", object@omega, "\n\n")
            cat("This means that at each point, the top", object@omega * 100,
                "percent of models in the last", object@lambda,
                "observations are selected and weighted for prediction")
          }
)


setMethod("show",
          signature("dets_hat"),
          function(object) {
            cat("## Predictions for DETS ##\n\n")
            cat("Check:\n")
            cat("slot @y_hat for the combined predictions\n")
            cat("slot @Y_hat for the predictions of each base model\n")
            cat("slot @Y_committee to see the selected models at each prediction\n")
            cat("slot @W for the weights of each base model\n")
          }
)


