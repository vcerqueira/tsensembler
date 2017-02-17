#' se
#'
#' Utility function to compute pointwise squared error (SE)
#' @param true True values
#' @param pred Predicted values
#' @export
se <- function(true, pred) {
  stopifnot(length(true) == length(pred),
            is.numeric(true),
            is.numeric(pred))

  (true - pred) ^ 2
}

#' mse
#'
#' Utility function to compute mean squared error (MSE)
#' @param true True values
#' @param pred Predicted values
#' @export
mse <- function(true, pred) mean(se(true, pred))

#' rmse
#'
#' Utility function to compute Root Mean Squared Error (RMSE)
#' @param true True values
#' @param pred Predicted values
#' @export
rmse <- function(true, pred) sqrt(mse(true, pred))

#' ae
#'
#' Absolute Error loss function. Element-wise computation
#' @param true True values
#' @param pred Predicted values
#' @export
ae <- function(true, pred) {
  stopifnot(length(true) == length(pred),
            is.numeric(true),
            is.numeric(pred))

  abs(true - pred)
}



#' mae
#'
#' Mean Absolute Error loss function.
#' @param true True values
#' @param pred Predicted values
#' @export
mae <- function(true, pred) mean(ae(true, pred))

#' sle
#'
#' Squared Logarithmic Error. Element-wise computation
#' @param true True values
#' @param pred Predicted values
#' @export
sle <- function(true, pred) {
  stopifnot(length(true) == length(pred),
            is.numeric(true),
            is.numeric(pred))

  (log(1 + true) - log(1 + pred)) ^ 2
}

#' msle
#'
#' Mean Squared Logarithmic Error. Element-wise computation
#' @param true True values
#' @param pred Predicted values
#' @export
msle <- function(true, pred) mean(sle(true, pred))

#' Log Percentual Difference
#'
#' This function computes the log percentual difference
#' between two numeric vector estimators
#'
#' @param x baseline estimator. Numeric vector
#' @param y Estimator to compare with the baseline. Numeric Vector
#'
#' @return Logarithmic percentual difference
#' between the two estimators
#'
#' @export
perf_diff <- function(x, y) {
  xy_diff <- ((y - x) / x)

  sign(xy_diff) * log(abs(xy_diff + 1))
}
