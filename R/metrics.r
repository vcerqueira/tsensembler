#' se
#'
#' Utility function to compute pointwise squared error (SE)
#'
#' @param y A numeric vector representing the actual values.
#' @param y_hat A numeric vector representing the forecasted values.
#'
#' @return squared error of forecasted values.
#'
#' @export
se <- function(y, y_hat) {
  stopifnot(length(y) == length(y_hat),
            is.numeric(y),
            is.numeric(y_hat))

  (y - y_hat) ^ 2
}

#' mse
#'
#' Utility function to compute mean squared error (MSE)
#'
#' @inheritParams se
#'
#' @export
mse <- function(y, y_hat) mean(se(y, y_hat))

#' rmse
#'
#' Utility function to compute Root Mean Squared Error (RMSE)
#'
#' @inheritParams se
#'
#' @export
rmse <- function(y, y_hat) sqrt(mse(y, y_hat))

#' ae
#'
#' Element-wise computation of the absolute error loss function.
#'
#' @inheritParams se
#'
#' @export
ae <- function(y, y_hat) {
  stopifnot(length(y) == length(y_hat),
            is.numeric(y),
            is.numeric(y_hat))

  abs(y - y_hat)
}

#' mae
#'
#' Mean Absolute Error loss function.
#'
#' @inheritParams se
#'
#' @export
mae <- function(y, y_hat) mean(ae(y, y_hat))

#' sle
#'
#' Element-wise computation of the squared logarithmic error
#'
#' @inheritParams se
#'
#' @export
sle <- function(y, y_hat) {
  stopifnot(length(y) == length(y_hat),
            is.numeric(y),
            is.numeric(y_hat))

  (log(1 + y) - log(1 + y_hat)) ^ 2
}

#' msle
#'
#' Mean squared logarithmic error.
#'
#' @inheritParams se
#'
#' @export
msle <- function(y, y_hat) mean(sle(y, y_hat))

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

#' R squared
#'
#' Also known as the coefficient of determination
#'
#' @inheritParams se
#'
#' @export
r_squared <- function(y, y_hat) {
  1 - (sum((y - y_hat )^2) / sum((y - mean(y))^2))
}

#' Symmetric Mean Absolute Percentage Error
#'
#' @inheritParams se
#'
#' @export
smape <- function(y, y_hat) {
  (100 / length(y)) * sum(abs(y_hat - y) / ((abs(y_hat) + abs(y)) / 2))
}

#' Mean absolute Scaled Error
#'
#' Mean absolute scaled error using the previous value as
#' baseline (i.e. a random walk)
#'
#' @inheritParams se
#'
#' @export
mase <- function(y, y_hat) {
  len <- length(y)
  baseline <- c(NA_real_, y[-length(y)])
  den <- (len / (len-1)) * sum(ae(y, baseline), na.rm = TRUE)
  sum(abs(y - y_hat)) / den
}
