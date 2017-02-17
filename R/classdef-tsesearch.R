setOldClass("xts")
setOldClass("zoo")
setClassUnion("tsSource", c("list", "xts", "zoo"))
setClassUnion("OptionalNumeric", c("numeric","NULL"))
setClassUnion("OptionalCharacter", c("character","NULL"))
setClassUnion("OptionalLogical", c("logical","NULL"))
setClassUnion("OptionalList", c("list","NULL"))


#' tsesearch
#'
#' An S4 class that combines the meta-data about the Time Series Ensemble.
#' It wraps the ensemble procedure on top of the performance estimation workflows.
#'
#' @slot timeseries A time series or list of time series of class \code{tsSource}.
#' This class inherits the classes \code{xts}, \code{zoo} and \code{list}.
#'
#' @slot workflow Object of class Workflow. Check \code{\link[performanceEstimation]{Workflow-class}} for
#' further details.
#'
#' @slot embedding.dimension The maximum embedding dimension used to transform the series.
#' If \code{timeseries} is a single time series and \code{embedding.dimension} is \emph{NULL}
#' the parameter will be estimated using function.
#'
#' @slot nReps Number of Monte Carlo repetitions to perform the experiment simulations.
#' Defaults to 10. More information in \code{\link[performanceEstimation]{MonteCarlo-class}}.
#'
#' @slot szTrain Ratio of the available window size used for training.
#' More information at \code{\link[performanceEstimation]{MonteCarlo-class}}
#'
#' @slot szTest Ratio of the available window size used for testing.
#' More information at \code{\link[performanceEstimation]{MonteCarlo-class}}
#'
#' @slot modelOutput Logical. If TRUE, the methods applied to \code{tsensembler-class} return
#' a predictive model ensemble. Otherwise, it returns performance statistics based on
#' \code{\link[performanceEstimation]{performanceEstimation}}.
#'
#' @seealso \code{\link[performanceEstimation]{performanceEstimation}} for the ensemble estimation
#' procedure and \code{\link[performanceEstimation]{MonteCarlo-class}} for Monte Carlo simulation setup.
#' The methods \code{\link{ensembler}} or \code{\link{Multi.ensembler}} are used to run an
#' object of class \code{tsensembler}. The normal output of the previously mentioned methods are
#' of class \code{\link{tseModel-class}}. 
#'
#' @examples
#' \dontrun{
#' TSEr <- tsensembler(timeseries = ts,
#'                     workflow = Workflow(wf = 'TSE',
#'                                         learner = c('MARS', 'PPR'),
#'                                         learner.pars = list(mars = list(nk = c(5, 2),
#'                                                                         degree= c(3, 4)),
#'                                                             ppr = list(nterms = c(2,3,4))),
#'                                         varying.embed = TRUE,
#'                                         varying.trainwindow = FALSE,
#'                                         committee.ratio = .1,
#'                                         combine.type = "regret_theorem23",
#'                                         verbose = FALSE),
#'                     modelOutput = TRUE,
#'                     embedding.dimension = 30,
#'                     nReps = 2)
#' }
#'
#' @export
setClass("tsesearch",
         slots = c(timeseries = "tsSource",
                   workflow = "list",
                   embedding.dimension = "OptionalNumeric",
                   modelOutput = "logical",
                   nReps = "numeric",
                   szTrain = "numeric",
                   szTest = "numeric")
)

#' tsesearch
#'
#' An S4 class that combines the meta-data about the Time Series Ensemble.
#' It wraps the ensemble procedure on top of the performance estimation workflows.
#'
#' @param timeseries A time series or list of time series of class \code{tsSource}.
#' This class inherits the classes \code{xts}, \code{zoo} and \code{list}.
#'
#' @param workflow Object of class Workflow. Check \code{\link[performanceEstimation]{Workflow-class}} for
#' further details.
#'
#' @param embedding.dimension The maximum embedding dimension used to transform the series.
#' If \code{timeseries} is a single time series and \code{embedding.dimension}
#'
#' @param nReps Number of Monte Carlo repetitions to perform the experiment simulations.
#' Defaults to 10. More information in \code{\link[performanceEstimation]{MonteCarlo-class}}.
#'
#' @param szTrain Ratio of the available window size used for training.
#' More information at \code{\link[performanceEstimation]{MonteCarlo-class}}
#'
#' @param szTest Ratio of the available window size used for testing.
#' More information at \code{\link[performanceEstimation]{MonteCarlo-class}}
#'
#' @param modelOutput Logical. If TRUE, the methods applied to \code{tsensembler-class} return
#' a predictive model ensemble. Otherwise, it returns performance statistics based on
#' \code{\link[performanceEstimation]{performanceEstimation}}.
#'
#' @seealso \code{\link[performanceEstimation]{performanceEstimation}} for the ensemble estimation
#' procedure and \code{\link[performanceEstimation]{MonteCarlo-class}} for Monte Carlo simulation setup.
#' The methods \code{\link{ensembler}} or \code{\link{Multi.ensembler}} are used to run an
#' object of class \code{tsensembler}. The normal output of the previously mentioned methods are
#' of class \code{\link{tseModel-class}}. 
#'
#' @examples
#' \dontrun{
#' TSEr <- tsensembler(timeseries = ts,
#'                     workflow = Workflow(wf = 'TSE',
#'                                         learner = c('MARS', 'PPR'),
#'                                         learner.pars = list(mars = list(nk = c(5, 2),
#'                                                                         degree= c(3, 4)),
#'                                                             ppr = list(nterms = c(2,3,4))),
#'                                         varying.embed = TRUE,
#'                                         varying.trainwindow = FALSE,
#'                                         committee.ratio = .1,
#'                                         combine.type = "regret_theorem23",
#'                                         verbose = FALSE),
#'                     modelOutput = TRUE,
#'                     embedding.dimension = 30,
#'                     nReps = 2)
#' }
#'
#' @export
tsesearch <- function(timeseries,
                      workflow,
                      embedding.dimension = NULL,
                      modelOutput = FALSE,
                      nReps = 7L,
                      szTrain = .6,
                      szTest = .25) {
  if (missing(timeseries) || missing(workflow)) {
    stop('\nYou need to provide a time series, a workflow and an embedding dimension.\n')
  }

  .learners <- c("baggedtrees", "SVM", "FFNN", "RandomForest",
                 "MARS", "Cubist", "PPR", "GLM", "GBM", "SAE", "GP")

  workflow <- c(workflow)

  if (is.null(embedding.dimension)) {
    if (!is.list(timeseries)) {
      embedding.dimension <- 15
      message("needs refactor")
    } else {
      stop("can't estimate embedding dimension for multiple time series.
           specify a numeric embedding dimension or estimate it using k_hat function.")
    }
    }

  for (i in seq_along(workflow)) {
    if (workflow[[i]]@func %in% "TSE") {
      pars <- workflow[[i]]@pars
      pars$modelOutput <- modelOutput
      workflow[[i]]@name <- paste0("WF", i)
      L <- pars$learner
      Lpars <- pars$learner.pars

      if (!is.null(L) && any(!L %in% .learners)) {
        stop("Invalid \"learner\" parameter.", call. = FALSE)
      }
      if (is.null(L)) stop("Please provide a learner for all workflows. Type help(ensembler) for details", .call = FALSE)

      if (length(L) == 1L && L == "baggedtrees") next

      if (is.null(Lpars)) message("Using default learner parameters.")
      if (!pars$aggregationFUN %in% c("emase-scommittee", "emase-wcommittee",
                                    "regret-scommittee", "regret-wcommittee",
                                    "emase-all", "regret-all",
                                    "static-s", "static-w", "npmase")) {
        stop("Model combination type (\"aggregationFUN\" parameter) is invalid.")
      }
      if (pars$aggregationFUN %in% c("static-s", "static-w")) next
      cr <- pars$committee.ratio
      if (cr <= 0 || cr > 1) {
        stop("Please choose a valid top committee ratio (double between 0 and 1).")
      }
    }
  }

  if (length(embedding.dimension) > 1L) {
    message("Estimating a model for each embedding dimension.")
  }

  new("tsesearch",
      timeseries = timeseries,
      workflow = workflow,
      embedding.dimension = embedding.dimension,
      modelOutput = modelOutput,
      nReps = nReps,
      szTrain = szTrain,
      szTest = szTest)
}
