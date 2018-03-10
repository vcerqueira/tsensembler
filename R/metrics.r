#' Computing the squared error
#'
#' Utility function to compute pointwise squared error (SE)
#'
#' @param y A numeric vector representing the actual values.
#' @param y_hat A numeric vector representing the forecasted values.
#'
#' @return squared error of forecasted values.
#'
#' @family error/performance functions
#'
#' @keywords internal
#'
#' @export
se <- function(y, y_hat, ...) {
  stopifnot(length(y) == length(y_hat),
            is.numeric(y),
            is.numeric(y_hat))

  (y - y_hat) ^ 2
}

#' Computing the mean squared error
#'
#' Utility function to compute mean squared error (MSE)
#'
#' @inheritParams se
#'
#' @family error/performance functions
#'
#' @keywords internal
#'
#' @export
mse <- function(y, y_hat) mean(se(y, y_hat), na.rm = TRUE)

#' Computing the root mean squared error
#'
#' Utility function to compute Root Mean Squared Error (RMSE)
#'
#' @inheritParams se
#'
#' @family error functions
#'
#' @keywords internal
#'
#' @export
rmse <- function(y, y_hat) sqrt(mse(y, y_hat))

#' Computing the absolute error
#'
#' Element-wise computation of the absolute error loss function.
#'
#' @inheritParams se
#'
#' @family error/performance functions
#'
#' @keywords internal
#'
#' @export
ae <- function(y, y_hat, ...) {
  stopifnot(length(y) == length(y_hat),
            is.numeric(y),
            is.numeric(y_hat))

  abs(y - y_hat)
}

#' Computing the mean absolute error
#'
#' Mean Absolute Error loss function.
#'
#' @inheritParams se
#'
#' @family error/performance functions
#'
#' @keywords internal
#'
#' @export
mae <- function(y, y_hat) mean(ae(y, y_hat), na.rm = TRUE)

#' Computing R squared
#'
#' Also known as the coefficient of determination
#'
#' @inheritParams se
#'
#' @family error/performance functions
#'
#' @keywords internal
#'
#' @export
r_squared <- function(y, y_hat) {
  1 - (sum((y - y_hat)^2, na.rm=TRUE) / sum((y - mean(y, na.rm=TRUE))^2, na.rm=TRUE))
}

#' Computing the mean absolute scaled error
#'
#' Mean absolute scaled error using the previous value as
#' baseline (i.e. a random walk)
#'
#' @inheritParams se
#'
#' @family error/performance functions
#'
#' @keywords internal
#'
#' @export
mase <- function(y, y_hat) {
  len <- length(y)
  baseline <- c(NA_real_, y[-length(y)])
  den <- (len / (len-1)) * sum(ae(y, baseline), na.rm = TRUE)
  sum(abs(y - y_hat)) / den
}
