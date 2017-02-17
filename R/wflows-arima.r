#' arimaWorkflow
#'
#' This function implements a workflow of autoregressive and moving average models.
#' It is based on \code{auto.arima} function of \code{forecast} package, which
#' performs an optimal parameter tuning.
#'
#' @param form Formula, e.g. \code{target ~.}
#' @param train embedded time series used for training the base learners
#' @param test embedded time series used for testing
#' @param modelOutput Logical. If TRUE, the methods applied to \code{tsensembler-class} return
#' a predictive model ensemble. Otherwise, it returns performance statistics based on
#' \code{\link[performanceEstimation]{performanceEstimation}}.
#' @param ... Further parameters to pass to \code{tse.arima}
#'
#' @seealso \code{\link[forecast]{auto.arima}} for an automatic arima procedure.
#'
#' @return a list containing the true values and predicted values for test set
#'
#' @importFrom forecast auto.arima
#' @importFrom forecast forecast.Arima
#' @export
tse.arima <- function(form, train, test, modelOutput = FALSE, ...) {
  unembedded.series <- unembed.timeseries(train, test)

  train <- unembedded.series$train
  test  <- trues <- unembedded.series$test
  data  <- unembedded.series$data

  seq. <- seq_along(test)
  trainset <- train
  y_hat <- vnapply(seq., function(j) {
    trainset <- c(trainset, test[seq_len(j) - 1])
    model <- auto.arima(trainset)
    point.forecast <- forecast.Arima(model, h = 1)$mean

    as.vector(point.forecast)
  })

  res <- list(trues = trues, preds = y_hat)
  if (modelOutput) res <- c(res, model = auto.arima(data))

  res
}

#' ARCH(1,1)-ARMA process workflow
#'
#' The forecast is calculated in a point-wise fashion. That is, the process
#' is estimated at each step.
#'
#' @param form Formula
#' @param train embedded time series used for training the base learners
#' @param test embedded time series used for testing
#' @param modelOutput Logical. If TRUE, the methods applied to \code{tsensembler-class} return
#' a predictive model ensemble. Otherwise, it returns performance statistics based on
#' \code{\link[performanceEstimation]{performanceEstimation}}.
#' @param ... Further parameters to \code{arch11armaWorkflow}
#'
#' @importFrom rugarch ugarchfit
#' @importFrom rugarch ugarchspec
#' @importFrom rugarch ugarchforecast
#' 
#' @export
tse.arch11arma <- function(form, train, test, modelOutput = FALSE, ...) {
  embedding.dimension <- NCOL(train)
  unembedded.series <- unembed.timeseries(train, test)
  train <- unembedded.series$train
  test  <- unembedded.series$test
  data  <- unembedded.series$data

  preds <- numeric(length(test))
  for (i in seq_along(preds)) {
    arimamodel <- forecast::auto.arima(train)
    ar <- arimamodel$arma[1]
    ma <- arimamodel$arma[2]

    archarma <- rugarch::ugarchfit(
      rugarch::ugarchspec(mean.model = list(armaOrder=c(ar, ma)),
                          distribution.model = "std"),
      train)

    preds[i] <- rugarch::ugarchforecast(archarma, n.ahead = 1)@forecast$seriesFor[[1]]
    train <- c(train, test[i])
  }

  res <- list(trues = test, preds = preds)
  if (modelOutput) res <- c(res, model = archarma)

  res
}

#' auto arima with external features
#'
#' @param form formula
#' @param train train set
#' @param test test set
#' @param xreg.id id of external feats
#' @param ... Further parameters to the workflow
#'
#' @importFrom forecast auto.arima
#' @importFrom forecast Arima
#'
#' @export
tse.arima_xreg <- function(form, train, test, xreg.id, ...) {
  tgt <- get_target(form)
  y_tr <- train[ ,tgt]
  Xext_tr <- as.matrix(train[, xreg.id])
  n.tr <- nrow(train)
  data <- rbind.data.frame(train, test)

  Y <- data[ ,tgt]
  X_ext <- as.matrix(data[ ,xreg.id])

  model <- auto.arima(y = y_tr, xreg = Xext_tr)

  predictions <- fitted(
    Arima(y = Y,
          xreg = X_ext,
          model = model)
  )[(n.tr + 1):NROW(data)]

  trues <- Y[(n.tr + 1):NROW(data)]

  res <- list(trues = trues, preds = predictions)

  res
}



