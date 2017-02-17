setMethod("show",
          signature("tseModel"),
          function(object) {
            cat("Time Series Ensemble Model\n")
            cat("Ensemble composed by", object@N, "base Models.\n")
            cat("The distribution of the base Models is the following:\n")
            print(object@modelDist)
            cat("The formula is:\n")
            print(object@form)
            cat("The Workflow setup:\n")
            cat("MASE n:", object@ma.N,"\n")
            cat("Embedding Dimension:", object@embedding.dimension, "\n")
            cat("Committee Ratio:", object@committee.ratio,"\n")
            cat("Combine Type:", object@aggregationFUN, "\n")
            cat("Workflow Name:", object@wfName,"\n")
          })

#' Predicting on new data with a tseModel class object
#'
#' This is a \code{predict} method for predicting new data points using a
#' \code{tseModel} class object - refering to an ensemble for time series
#' forecasting.
#'
#' @seealso \code{\link{tseModel-class}} for details about the parameter setting
#' of the \code{object}, \code{\link{computePredictions}} for details on
#' the computation of predictions of the base learning models;
#' \code{\link{LearnerAnalysis}} for details on the weighting of models along the
#' test observations; \code{\link{combinePredictions}} on details on how the
#' base learners' predictions are combined.
#'
#' @param object A \strong{tseModel-class} object.
#' @param newdata An embedded time series with new data points to predict.
#'
#' @return predict method for \code{tseModel} class returns four objects:
#' \describe{
#'    \item{base.predictions}{A list containing the predictions of the base learners}
#'    \item{predictions}{The final predictions}
#'    \item{baseScore}{The weights of the base learners along the new data}
#'    \item{topBaseLearners}{The top \code{committee.ratio} learners along the new data}
#' }
#'
#' @examples
#'
#' \dontrun{
#' require(xts)
#' ts <- as.xts(rnorm(100L), order.by = Sys.Date() + rnorm(100L))
#' TSEr <- tsensembler(timeseries = ts,
#'                     form = target ~.,
#'                     workflow = Workflow(wf = 'TSE',
#'                                         learner = c('MARS', 'PPR'),
#'                                         learner.pars = list(mars = list(nk = c(5, 2),
#'                                                                         degree= c(3, 4)),
#'                                                             ppr = list(nterms = c(2,3,4))),
#'                                         varying.embed = TRUE,
#'                                         varying.trainwindow = FALSE,
#'                                         committee.ratio = .1,
#'                                         aggregationFUN = "emase-scommittee",
#'                                         verbose = FALSE),
#'                     embedding.dimension = 30,
#'                     nReps = 2,
#'                     modelOutput = TRUE)
#' tse.model <- ensembler(TSEr)
#' newdata <- embedStats(embed.timeseries(ts, embedding.dimension))
#'
#' preds <- predict(tse.model, newdata)
#' }
#'
#' @export
setMethod("predict",
          signature("tseModel"),
          function(object, newdata) {
            .wf <- object@wfName

            if (.wf == "SHoW") {
              nullcheck <- vlapply(object@baseModels, is.null)
              if (any(nullcheck)) {
                null.models <- which(nullcheck)
                models[null.models] <- NULL
              }
              predictions <- sapply(object@baseModels, function(model) predict(model, newdata))

              if (is.null(dim(predictions))) {
                preds <- mean(predictions)
              } else {
                preds <- rowMeans(predictions)
              }
              preds <- tse_hat(y_hat = preds)
            } else {
              Scoring <- NULL
              seq. <- seq_len(object@N)
              target <- get_target(object@form)

              predictions <- unlist(lapply(seq., function(.i) {
                computePredictions(.Model = object@baseModels[.i],
                                   data = newdata,
                                   form = object@form)
              }), recursive = FALSE)

              if (.wf %in% "DHeW") {
                if (!is.numeric(object@ma.N)) stop("ma.N must be numeric.")

                Scoring <- LearnerAnalysis(predictions = predictions,
                                           ma.N = object@ma.N,
                                           committee.ratio = object@committee.ratio,
                                           preweights = object@preWeights,
                                           aggregationFUN = object@aggregationFUN,
                                           Mstar = object@Mstar)
              }
              preds <- combinePredictions(predictions = predictions,
                                          test = newdata,
                                          topmodels = Scoring$Toplearners,
                                          modelscores = Scoring$fScores,
                                          modelpreweights = object@preWeights,
                                          aggregationFUN = object@aggregationFUN)

              preds <- tse_hat(y_hat = preds, 
                               Y_hat = predictions,
                               Y_W = Scoring$fScores,
                               Y_committee = Scoring$Toplearners,
                               Y_rank = NULL,
                               data = newdata)
            }
            preds
          })


#' Forecasting with a Dynamic Ensemble for Time Series tasks.
#'
#' \code{forecast} is a method under \code{tseModel-class} domain that
#' predicts the next \code{h} data points given a \code{tseModel-class}
#' trained on a given \code{timeseries}.
#'
#' @seealso \code{\link{predict}} for a closely related method.
#'
#' @param object \code{tseModel-class} object.
#' @param timeseries Time series used to train the \code{object}.
#' @param h Number of periods for forecasting.
#'
#' @return \code{h} forecasting points ahead of a given \code{timeseries}.
#'
#' @examples
#'
#' \dontrun{
#' require(xts)
#' ts <- as.xts(rnorm(100L), order.by = Sys.Date() + rnorm(100L))
#' TSEr <- tsensembler(timeseries = ts,
#'                     form = target ~.,
#'                     workflow = Workflow(wf = 'TSE',
#'                                         learner = c('MARS', 'PPR'),
#'                                         learner.pars = list(mars = list(nk = c(5, 2),
#'                                                                         degree= c(3, 4)),
#'                                                             ppr = list(nterms = c(2,3,4))),
#'                                         varying.embed = TRUE,
#'                                         varying.trainwindow = FALSE,
#'                                         committee.ratio = .1,
#'                                         aggregationFUN = "emase-wcommittee",
#'                                         verbose = FALSE),
#'                     embedding.dimension = 30,
#'                     nReps = 2,
#'                     modelOutput = TRUE)
#' tse.model <- ensembler(TSEr)
#'
#' pointforecast <- forecast(tse.model, ts, h = 4)
#' }
#'
#' @export
setGeneric("forecast", function(object, timeseries, h) {
  standardGeneric("forecast")
})

#' Forecasting with a Dynamic Ensemble for Time Series tasks.
#'
#' \code{forecast} is a method under \code{tseModel-class} domain that
#' predicts the next \code{h} data points given a \code{tseModel-class}
#' trained on a given \code{timeseries}.
#'
#' @seealso \code{\link{predict}} for a closely related method.
#'
#' @param object \code{tseModel-class} object.
#' @param timeseries Time series used to train the \code{object}.
#' @param h Number of periods for forecasting.
#'
#' @return \code{h} forecasting points ahead of a given \code{timeseries}.
#'
#' @examples
#'
#' \dontrun{
#' require(xts)
#' ts <- as.xts(rnorm(100L), order.by = Sys.Date() + rnorm(100L))
#' TSEr <- tsensembler(timeseries = ts,
#'                     form = target ~.,
#'                     workflow = Workflow(wf = 'TSE',
#'                                         learner = c('MARS', 'PPR'),
#'                                         learner.pars = list(mars = list(nk = c(5, 2),
#'                                                                         degree= c(3, 4)),
#'                                                             ppr = list(nterms = c(2,3,4))),
#'                                         varying.embed = TRUE,
#'                                         varying.trainwindow = FALSE,
#'                                         committee.ratio = .1,
#'                                         aggregationFUN = "emase-scommittee",
#'                                         verbose = FALSE),
#'                     embedding.dimension = 30,
#'                     nReps = 2,
#'                     modelOutput = TRUE)
#' tse.model <- ensembler(TSEr)
#'
#' pointforecast <- forecast(tse.model, ts, h = 4)
#' }
#'
#' @export
setMethod("forecast",
          signature("tseModel"),
          function(object, timeseries, h) {
            .wf <- object@wfName

            K <- object@embedding.dimension
            len <- NROW(timeseries)

            rawsignal <- as.vector(timeseries)[len:(len-(K-2L))]
            forecasts <- numeric(h)
            for (t in seq_len(h)) {
              signal <- readSignal(rawsignal)
              pointforecast <- predict(object, signal)$predictions
              forecasts[t]  <- pointforecast
              rawsignal    <- FIFO(rawsignal, pointforecast)
            }
            names(forecasts) <- paste("t", seq_along(forecasts), sep = "+")
            forecasts
          })
