setOldClass("xts")
setOldClass("zoo")
setClassUnion("embeddedseriesSource", c("data.frame", "xts", "zoo"))
setClassUnion("OptionalNumeric", c("numeric","NULL"))
setClassUnion("OptionalCharacter", c("character","NULL"))
setClassUnion("OptionalLogical", c("logical","NULL"))
setClassUnion("OptionalList", c("list","NULL"))


#' Time Series Ensemble
#'
#' \code{tsensembler} identifies the model specs for creating the dynamic
#' heterogeneous ensemble for time series forecasting.
#'
#' @slot timeseries A time series of class xts
#' @slot learner base models to learn the data
#' @slot learner.pars list with the parameter setting of \code{learner}
#' @slot varying.embed Logical. If TRUE, each learning model in \code{learner}
#' is trained using three different subsets of the data w.r.t. embedding dimension.
#' @slot varying.trainwindow Logical. If TRUE, each learning model in \code{learner}
#' is trained using three different subsets of the data w.r.t. training window.
#' @slot committee.ratio Double between 0 and 1 representing the ratio of
#' learners that should be considered at each prediction time.
#' @slot embedding.dimension The maximum embedding dimension used to transform the series.
#' If \code{timeseries} is a single time series and \code{embedding.dimension} is \emph{NULL}
#' the parameter will be estimated using false nearest neighbors.
#' @slot ma.N The number of periods to average over the Squared Error when
#' computing \strong{MASE} (\emph{Moving Average Squared Error}). This parameter is only
#' used when the function applied to weight the base models is based on
#' MASE. The rationale behind this simple heuristic is that base models are weighted according
#' to their recent performance. Here, recent performance is formalized and quantified as
#' inversely proportional to the models' MASE. The dynamics of the moving average
#' yield a flexibility to the combined model, in the sense that it is self-adaptable when
#' concept drift occurs. Particularly, this parameter \code{ma.N} controls
#' the reactiveness of the learning system to such events.
#' A small value of P leads to greater reactiveness, but also makes the ensemble
#' susceptible to be fooled by outliers. Conversely, higher values of \code{ma.N}
#' lead to greater stability, while sacrificing some responsiveness. This trade-off
#' is known in the literature as the stability-plasticity dilemma (Carpenter et al., 1991).
#' @slot aggregationFUN The function name used to combine the base learners. See \link{combinePredictions}
#' for a comprehensive explanation.
#' @slot verbose Logical. If TRUE, information about the learning procedure status is printed into the console.
#'
#' @seealso \code{\link{tsesearch}}
#' @export
setClass("tsensembler",
         slots = c(timeseries = "embeddedseriesSource",
                   learner = "OptionalCharacter",
                   learner.pars = "OptionalList",
                   varying.embed = "OptionalLogical",
                   varying.trainwindow = "OptionalLogical",
                   committee.ratio = "OptionalNumeric",
                   embedding.dimension = "OptionalNumeric",
                   ma.N = "OptionalNumeric",
                   aggregationFUN = "OptionalCharacter",
                   verbose = "logical")
)

#' Time Series Ensemble
#'
#' \code{tsensembler} identifies the model specs for creating the dynamic
#' heterogeneous ensemble for time series forecasting.
#'
#' @param timeseries A time series of class xts
#' @param learner base models to learn the data
#' @param learner.pars list with the parameter setting of \code{learner}
#' @param varying.embed Logical. If TRUE, each learning model in \code{learner}
#' is trained using three different subsets of the data w.r.t. embedding dimension.
#' @param varying.trainwindow Logical. If TRUE, each learning model in \code{learner}
#' is trained using three different subsets of the data w.r.t. training window.
#' @param committee.ratio Double between 0 and 1 representing the ratio of
#' learners that should be considered at each prediction time.
#' @param embedding.dimension The maximum embedding dimension used to transform the series.
#' If \code{timeseries} is a single time series and \code{embedding.dimension} is \emph{NULL}
#' the parameter will be estimated using false nearest neighbors method.
#' @param ma.N The number of periods to average over the Squared Error when
#' computing \strong{MASE} (\emph{Moving Average Squared Error}). This parameter is only
#' used when the function applied to weight the base models is based on
#' MASE. The rationale behind this simple heuristic is that base models are weighted according
#' to their recent performance. Here, recent performance is formalized and quantified as
#' inversely proportional to the models' MASE. The dynamics of the moving average
#' yield a flexibility to the combined model, in the sense that it is self-adaptable when
#' concept drift occurs. Particularly, this parameter \code{ma.N} controls
#' the reactiveness of the learning system to such events.
#' A small value of P leads to greater reactiveness, but also makes the ensemble
#' susceptible to be fooled by outliers. Conversely, higher values of \code{ma.N}
#' lead to greater stability, while sacrificing some responsiveness. This trade-off
#' is known in the literature as the stability-plasticity dilemma (Carpenter et al., 1991).
#' @param aggregationFUN The function name used to combine the base learners. See \link{combinePredictions}
#' for a comprehensive explanation.
#' @param verbose Logical. If TRUE, information about the learning procedure status is printed into the console.
#'
#' @seealso \code{\link{tsesearch}}
#' @export
tsensembler <- function(timeseries,
                        learner = NULL,
                        learner.pars = NULL,
                        varying.embed = FALSE,
                        varying.trainwindow = FALSE,
                        committee.ratio = .1,
                        embedding.dimension = 40,
                        ma.N = 50,
                        aggregationFUN = "emase-wcommittee",
                        verbose = FALSE) {
  if (missing(timeseries)) {
    stop('\nYou need to provide a time series of class \'xts\'.\n')
  }

  .learners <- c("baggedtrees", "SVM", "FFNN", "RandomForest",
                 "MARS", "Cubist", "PPR", "GLM", "GBM", "SAE", "GP")

  if (!any(aggregationFUN %in% c("emase-scommittee", "emase-wcommittee",
                                 "regret-scommittee", "regret-wcommittee",
                                 "emase-all", "regret-all",
                                 "static-s", "static-w", "npmase"))) {
    stop("Model combination type (\"aggregationFUN\" parameter) is invalid.")
  }

  if (!is.null(learner) && any(!learner %in% .learners)) {
    stop("Invalid \"learner\" parameter. Check list of available learners.")
  }

  if (is.null(learner.pars)) cat("Using default learner parameters.\n")

  if (is.null(embedding.dimension) && !is.data.frame(timeseries)) {
    cat("Estimating Embedding Dimension by Auto-correlation function.\n")
    embedding.dimension <- 15
    message("needs refactor")
  }

  if (any(committee.ratio <= 0) || any(committee.ratio > 1)) {
    stop("Please choose a valid top committee ratio (double between 0 and 1).")
  }

  new("tsensembler",
      timeseries = timeseries,
      learner = learner,
      learner.pars = learner.pars,
      varying.embed = varying.embed,
      varying.trainwindow = varying.trainwindow,
      committee.ratio = committee.ratio,
      embedding.dimension = embedding.dimension,
      ma.N = ma.N,
      aggregationFUN = aggregationFUN,
      verbose = verbose)
}
