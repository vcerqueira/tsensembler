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

#' @rdname forecast
setMethod("forecast",
          signature("DETS"),
          function(object, h) {

            forecasts <- numeric(h)
            for (o in seq_len(h)) {
              p_forecast <- forecast_dets(object)
              forecasts[o] <- p_forecast

              object@recent_series <-
                rbind.data.frame(object@recent_series,
                                 struc_embed(object, p_forecast))
            }
            names(forecasts) <- paste("t", seq_along(forecasts), sep = "+")

            forecasts
          })

forecast_dets <-
  function(object) {
    next_embed <- struc_embed(object, -1)

    Y_hat <- predict(object@base_ensemble, next_embed)

    Y_hat_recent <-
      predict(object@base_ensemble, object@recent_series)
    Y_recent <- get_y(object@recent_series, object@form)

    scores <-
      model_recent_performance(Y_hat_recent,
                               Y_recent,
                               object@lambda,
                               object@omega,
                               object@base_ensemble@pre_weights)

    model_scores <- scores$model_scores
    top_models <- scores$top_models

    W <- as.vector(model_scores[nrow(model_scores), ])
    C <- top_models[[length(top_models)]]

    Y_hat_j <- Y_hat[, C]
    W_j <- proportion(W[C])

    point_forecast <- sum(Y_hat_j * W_j)

    point_forecast
  }
