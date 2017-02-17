setMethod("show",
          signature("tsensembler"),
          function(object) {
            series <- deparse(substitute(object@timeseries))

            cat("tsensembler class object \n",
                "Use ensembler-method to train a dynamic ensemble with this
                object into", series, "timeseries\n\n")
            cat("Experimental Setup for Model Estimation:\n")
            cat("Learner:", object@learner, "\n")
            cat("Embedding dimension:", object@embedding.dimension, "\n")
            cat("Varying Embed:", object@varying.embed, "\n")
            cat("Varying Train Window:", object@varying.trainwindow, "\n")
            cat("Committee Ratio:", object@committee.ratio, "\n")
            cat("Aggregation FUN:", object@aggregationFUN, "\n")
            cat("Moving Average Periods:", object@ma.N, "\n\n")

            cat("Available methods:", ifelse(is.data.frame(object@timeseries),
                                             "ensembler",
                                             "univariate.ensembler"), "\n")
          })

#' Learning a Dynamic Ensemble for Time Series Forecasting
#'
#' This is the implementation of the learning procedure of the ensemble method.
#' The application is run on an object of class \code{tsensembler}, which provide
#' the options about how the ensemble should be construct such as: the embedding
#' dimension; the number of Monte Carlo repetitions, etc.
#' This method differs from \code{gs.ensembler} in the fact that \emph{ensembler}
#' does not estimate several models. \emph{ensembler} simply trains an Dynamic
#' Ensemble according to the specs in the object \strong{obj} of class \code{tsensembler}.
#'
#' @param obj An object of class \code{tsensembler}.
#'
#' @return A predictive model ensemble of class \code{tseModel}.
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
setGeneric("univariate.ensembler", function(obj) {
  standardGeneric("univariate.ensembler")
})

#' Learning a Dynamic Ensemble for Time Series Forecasting
#'
#' This is the implementation of the learning procedure of the ensemble method.
#' The application is run on an object of class \code{tsensembler}, which provide
#' the options about how the ensemble should be construct such as: the embedding
#' dimension; the number of Monte Carlo repetitions, etc.
#' This method differs from \code{gs.ensembler} in the fact that \emph{ensembler}
#' does not estimate several models. \emph{ensembler} simply trains an Dynamic
#' Ensemble according to the specs in the object \strong{obj} of class \code{tsensembler}.
#'
#' @param obj An object of class \code{tsensembler}.
#'
#' @return A predictive model ensemble of class \code{tseModel}.
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
setMethod("univariate.ensembler",
          signature("tsensembler"),
          function(obj) {
            embedded.series <- embed.timeseries(obj@timeseries, obj@embedding.dimension)
            seriestarget <- embedded.series[ , "target"]
            K <- obj@embedding.dimension
            row.size   <- NROW(embedded.series)
            dataStats <- embedStats(embedded.series[ ,!colnames(embedded.series) %in% "target"])
            eseries <- cbind.data.frame(target = seriestarget, dataStats)
            nstats <- NCOL(eseries) - K

            if (obj@varying.embed) {
              Ksplits <- c(K, K/2, K/4)
            } else {
              Ksplits <- K
            }
            if (obj@varying.trainwindow) {
              Rsplits <- floor(c(row.size, row.size/2, row.size/4)) - 1
            } else {
              Rsplits <- row.size
            }

            Learningensemble <- Learntse(eseries,
                                         target ~.,
                                         obj@learner,
                                         obj@learner.pars,
                                         Ksplits,
                                         Rsplits,
                                         nstats,
                                         K,
                                         obj@verbose)

            .tse <- tseModel(baseModels = Learningensemble$learningModels,
                             preWeights = Learningensemble$preweights,
                             form = target ~.,
                             ma.N = obj@ma.N,
                             embedding.dimension = K,
                             committee.ratio = obj@committee.ratio,
                             aggregationFUN = obj@aggregationFUN,
                             Mstar = Learningensemble$Mstar)
            .tse
          }
)

#' Generic Dynamic Ensembles for Parametric Learning Models
#'
#' ensembler is a method from class \code{tsensembler}. This method
#' applies learns a dynamic ensemble from the dataset provided.
#'
#' @param obj object from class \strong{tsensembler}
#' @param form formula
setGeneric("ensembler", function(obj, form) {
  standardGeneric("ensembler")
})

#' Generic Dynamic Ensembles for Parametric Learning Models
#'
#' ensembler is a method from class \code{tsensembler}. This method
#' applies learns a dynamic ensemble from the dataset provided.
#'
#' @param obj object from class \strong{tsensembler}
#' @param form formula
setMethod("ensembler",
          signature("tsensembler"),
          function(obj, form) {
            target <- .splitBy(deparse(form), " ")[1]
            trues  <- obj@timeseries[ , target]

            row.size   <- NROW(obj@timeseries)

            if (obj@varying.trainwindow) {
              Rsplits <- floor(c(row.size, row.size/2, row.size/4)) - 1
            } else {
              Rsplits <- row.size
            }

            Learningensemble <- Learntse(obj@timeseries,
                                         target ~.,
                                         obj@learner,
                                         obj@learner.pars,
                                         obj@embedding.dimension,
                                         Rsplits,
                                         0L,
                                         obj@embedding.dimension,
                                         obj@verbose)

            .tse <- tseModel(baseModels = Learningensemble$learningModels,
                             preWeights = Learningensemble$preweights,
                             form = target ~ .,
                             ma.N = obj@ma.N,
                             embedding.dimension = obj@embedding.dimension,
                             committee.ratio = obj@committee.ratio,
                             aggregationFUN = obj@aggregationFUN,
                             Mstar = Learningensemble$Mstar)
            .tse
          }
)
