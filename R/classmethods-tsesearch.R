#' Learning a Dynamic Ensemble for Time Series Forecasting
#'
#' This is the implementation of the learning procedure of the ensemble method.
#' The application is run on an object of class \code{tsensembler}, which provide
#' the options about how the ensemble should be construct such as: the embedding
#' dimension; the number of Monte Carlo repetitions, etc.
#'
#' This method is built on top of \code{\link[performanceEstimation]{performanceEstimation}}.
#' It provides a framework for estimating the performance of many workflows.
#'
#' @param obj An object of class \code{tsensembler}.
#'
#' @return If the \emph{modelOutput} parameter from \code{tsensembler-class} object is set to TRUE,
#' the function will return a predictive model ensemble of class \code{tseModel}. Otherwise,
#' performance metrics are return.
#'
#' @seealso \code{\link{tsensembler-class}} for details about the \strong{obj} and
#' \code{\link{tseModel-class}} for details on the generated model.
#' \code{\link[performanceEstimation]{MonteCarlo-class}} comprises an explanation
#' about the Monte Carlo simulation procedure. \code{\link{TSE}} for details on the workflow
#' for the time series ensemble.
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
#' }
#'
#' @export
setGeneric("gs.ensembler", function(obj) {
  standardGeneric("gs.ensembler")
})

#' Learning a Dynamic Ensemble for Time Series Forecasting
#'
#' This is the implementation of the learning procedure of the ensemble method.
#' The application is run on an object of class \code{tsensembler}, which provide
#' the options about how the ensemble should be construct such as: the embedding
#' dimension; the number of Monte Carlo repetitions, etc.
#'
#' This method is built on top of \code{\link[performanceEstimation]{performanceEstimation}}.
#' It provides a framework for estimating the performance of many workflows.
#'
#' @param obj An object of class \code{tsensembler}.
#'
#' @return If the \emph{modelOutput} parameter from \code{tsensembler-class} object is set to TRUE,
#' the function will return a predictive model ensemble of class \code{tseModel}. Otherwise,
#' performance metrics are return.
#'
#' @seealso \code{\link{tsensembler-class}} for details about the \strong{obj} and
#' \code{\link{tseModel-class}} for details on the generated model.
#' \code{\link[performanceEstimation]{MonteCarlo-class}} comprises an explanation
#' about the Monte Carlo simulation procedure. \code{\link{TSE}} for details on the workflow
#' for the time series ensemble.
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
#' }
#'
#' @export
setMethod("gs.ensembler",
          signature("tsesearch"),
          function(obj) {
            results <- list()
            for (K in obj@embedding.dimension) {
              embedded.series <- embed.timeseries(obj@timeseries, K)
              for (wf in obj@workflow) {
                wf@pars$embedding.dimension <- K
                estimateTSE <- performanceEstimation::performanceEstimation(
                  tasks     = PredTask(form = target ~.,
                                       data = embedded.series,
                                       copy = TRUE,
                                       type = "regr",
                                       taskName = paste0("task", "_K=", K)),
                  workflows = wf,
                  estTask   = EstimationTask(metrics = c("mse", "rmse", "mae"),
                                             method = MonteCarlo(obj@nReps, obj@szTrain, obj@szTest))
                )
                results <- c(results, list(estimateTSE))
              }
            }
            results <- .assemble(results)
            if (obj@modelOutput) results <- getTSE(results)
            results
          }
)


#' Learning a Dynamic Ensemble for Time Series Forecasting for \emph{Multi}ple Series
#'
#' This is the implementation of the learning procedure of the ensemble method for
#' \strong{Multi}ple Time Series.
#' The application is run on an object of class \code{tsensembler}, which provide
#' the options about how the ensemble should be construct such as: the embedding
#' dimension; the number of Monte Carlo repetitions, etc.
#'
#' This method is built on top of \code{\link[performanceEstimation]{performanceEstimation}}.
#' It provides a framework for estimating the performance of many workflows.
#'
#' @param obj An object of class \code{tsensembler}. The main difference between the
#' \code{ensembler} and \code{Multi.ensembler} methods with respect to \strong{obj} is that
#' the latter \emph{timeseries} parameter is a list of \code{xts} timeseries.
#' @param save.each Logical: Should each time series' experiment be saved? Defaults to TRUE
#' @param save.signature The name of the file template to save the experiments
#' @param auto.tune Logical: Should auto.ensemble be executed? Defaults to FALSE.
#' @param ... Further parameters to \code{auto.ensembler}
#'
#' @return If the \emph{modelOutput} parameter from \code{tsensembler-class} object is set to TRUE,
#' the function will return a predictive model ensemble of class \code{tseModel}. Otherwise,
#' performance metrics are return.
#'
#' @seealso \code{\link{tsensembler-class}} for details about the \strong{obj} and
#' \code{\link{tseModel-class}} for details on the generated model.
#' \code{\link[performanceEstimation]{MonteCarlo-class}} comprises an explanation
#' about the Monte Carlo simulation procedure. \code{\link{TSE}} for details on the workflow
#' for the time series ensemble.
#'
#' @examples
#'
#' \dontrun{
#' require(xts)
#' ts1 <- as.xts(rnorm(100L), order.by = Sys.Date() + rnorm(100L))
#' ts2 <- as.xts(rnorm(200L), order.by = Sys.Date() + rnorm(200L))
#' TSEr <- tsensembler(timeseries = list(one.ts = ts1, two.ts = ts2),
#'                     workflow = Workflow(wf = 'TSE',
#'                                         learner = c('MARS', 'PPR'),
#'                                         learner.pars = list(mars = list(nk = c(5, 2),
#'                                                                         degree= c(3, 4)),
#'                                                             ppr = list(nterms = c(2,3,4))),
#'                                         varying.embed = TRUE,
#'                                         varying.trainwindow = FALSE,
#'                                         committee.ratio = .1,
#'                                         aggregationFUN = "regret-wcommittee",
#'                                         verbose = FALSE),
#'                     embedding.dimension = 30,
#'                     nReps = 2,
#'                     modelOutput = TRUE)
#' tse.model <- ensembler(TSEr)
#' }
#'
#' @export
setGeneric("Multi.ensembler", function(obj,
                                       save.each = FALSE,
                                       save.signature = "Exp",
                                       auto.tune = FALSE, ...) {
  standardGeneric("Multi.ensembler")
})

#' Learning a Dynamic Ensemble for Time Series Forecasting for \emph{Multi}ple Series
#'
#' This is the implementation of the learning procedure of the ensemble method for
#' \strong{Multi}ple Time Series.
#' The application is run on an object of class \code{tsensembler}, which provide
#' the options about how the ensemble should be construct such as: the embedding
#' dimension; the number of Monte Carlo repetitions, etc.
#'
#' This method is built on top of \code{\link[performanceEstimation]{performanceEstimation}}.
#' It provides a framework for estimating the performance of many workflows.
#'
#' @param obj An object of class \code{tsensembler}. The main difference between the
#' \code{ensembler} and \code{Multi.ensembler} methods with respect to \strong{obj} is that
#' the latter \emph{timeseries} parameter is a list of \code{xts} timeseries.
#' @param save.each Logical: Should each time series' experiment be saved? Defaults to TRUE
#' @param save.signature The name of the file template to save the experiments
#' @param auto.tune Logical: Should auto.ensemble be executed? Defaults to FALSE.
#' @param ... Further parameters to \code{auto.ensembler}
#'
#' @return If the \emph{modelOutput} parameter from \code{tsensembler-class} object is set to TRUE,
#' the function will return a predictive model ensemble of class \code{tseModel}. Otherwise,
#' performance metrics are return.
#'
#' @seealso \code{\link{tsensembler-class}} for details about the \strong{obj} and
#' \code{\link{tseModel-class}} for details on the generated model.
#' \code{\link[performanceEstimation]{MonteCarlo-class}} comprises an explanation
#' about the Monte Carlo simulation procedure. \code{\link{TSE}} for details on the workflow
#' for the time series ensemble.
#'
#' @examples
#'
#' \dontrun{
#' require(xts)
#' ts1 <- as.xts(rnorm(100L), order.by = Sys.Date() + seq_len(100L))
#' ts2 <- as.xts(rnorm(200L), order.by = Sys.Date() + seq_len(200L))
#' TSEr <- tsensembler(timeseries = list(one.ts = ts1, two.ts = ts2),
#'                     workflow = Workflow(wf = 'TSE',
#'                                         learner = c('MARS', 'PPR'),
#'                                         learner.pars = list(mars = list(nk = c(5, 2),
#'                                                                         degree= c(3, 4)),
#'                                                             ppr = list(nterms = c(2,3,4))),
#'                                         varying.embed = TRUE,
#'                                         varying.trainwindow = FALSE,
#'                                         committee.ratio = .1,
#'                                         aggregationFUN = "regret-wcommittee",
#'                                         verbose = FALSE),
#'                     embedding.dimension = 30,
#'                     nReps = 2,
#'                     modelOutput = TRUE)
#' tse.model <- ensembler(TSEr)
#' }
#'
#' @export
setMethod("Multi.ensembler",
          signature("tsesearch"),
          function(obj,
                   save.each = FALSE,
                   save.signature = "Exp",
                   auto.tune = FALSE, ...) {
            if (!methods::is(obj@timeseries, "list")) {
              stop("Provide the time series as a list.\n")
            }
            seq. <- seq_along(obj@timeseries)
            if (is.null(names(obj@timeseries))) names(obj@timeseries) <- paste0("TS", seq.)
            ensembles <- list()
            for (ts in seq.) {
              cat(paste("Experimenting on TS", ts, "of", length(obj@timeseries), "\n"))
              tsensemblerObj <-   tsesearch(timeseries = obj@timeseries[[ts]],
                                            workflow = obj@workflow,
                                            embedding.dimension = obj@embedding.dimension,
                                            nReps = obj@nReps,
                                            szTrain = obj@szTrain,
                                            szTest = obj@szTest,
                                            modelOutput = obj@modelOutput)

              if (!auto.tune) {
                ensembles[[ts]] <- gs.ensembler(tsensemblerObj)
              } else {
                ensembles[[ts]] <- auto.ensembler(tsensemblerObj, ...)
              }
              if (save.each) {
                ExpRes <- ensembles[[ts]]
                save(ExpRes, file = paste0(save.signature, "_", names(obj@timeseries[ts]), ".Rdata"))
              }
            }
            ensembles
          }
)

setMethod("show",
          signature("tsesearch"),
          function(object) {
            cat("tsesearch class object \n\n")
            if (length(object@workflow) == 1L) {
              print(object@workflow[[1]])
            } else {
              print(object@workflow)
            }
            cat("Experimental Setup for Model Estimation:\n")
            cat("Embedding dimension:", object@embedding.dimension, "\n")
            cat("Number of Monte Carlo Repetitions:", object@nReps, "\n")
            cat("Train window size:", object@szTrain, "\n")
            cat("Test window size:", object@szTest, "\n\n\n")
            cat("Task is set for:",
                ifelse(object@modelOutput, "Modelling", "Estimation"), "\n\n")
          })
