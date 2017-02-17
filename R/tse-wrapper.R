#' Auto Tuned Time Series Ensemble
#'
#' This function is built on top of tsensembler. It is a \"for dummies\"
#' implementation that performs a grid search to get the best possible
#' model for the given time series.
#'
#' @param timeseries A time series of class xts.
#' @param search.level The search for the best model can be done in three levels:
#' \(i\) \emph{simple}; \(ii\) moderate; and \(iii\) extensive.
#' Beware that an \"extensive\" search will take a lot of time to run.
#'
#' @return The output is the best performing workflow as a predictive model.
#'
#' @seealso \code{\link{tseModel-class}} for details of the object class
#' returned. \code{\link[performanceEstimation]{performanceEstimation}} on details
#' of the implementation of the estimation procedure;
#' \code{\link[performanceEstimation]{MonteCarlo-class}} on details of the
#' Monte Carlo simulation method. \code{\link{TSE}} for details on the workflow
#' for the time series ensemble.
#'
#' @examples
#'
#' \dontrun{
#' require(xts)
#' ts <- as.xts(rnorm(100L), order.by = Sys.Date() + rnorm(100L))
#' tse <- auto.ensembler(ts, search.level = "simple")
#' }
#'
#' @export
#' @import performanceEstimation
auto.ensembler <- function(timeseries, search.level = "moderate") {
  if (!"xts" %in% class(timeseries)) {
    stop("Please provide a timeseries of \"xts\" class")
  }
  if (!search.level %in% c("simple", "moderate", "extensive")) {
    stop('\nProvide a valid search plan. Options are \"simple\",
         \"moderate\" and \"extensive\".\n')
  }

  series <- deparse(substitute(timeseries))

  if (search.level == "extensive") {
    tsObj <- tsesearch(timeseries = timeseries,
                       workflow =  c(Workflow(wf = 'arimaWorkflow'),
                                     Workflow(wf = 'TSE',
                                              learner <- c('MARS', 'GLM','FFNN',
                                                           'RandomForest', 'PPR',
                                                           'SVM', 'Cubist', 'GBM', 'GP'),
                                              learner.pars <- list(ffnn = list(size = c(10, 30, 50),
                                                                               decay = c(0, 0.05),
                                                                               maxit = 750),
                                                                   mars = list(nk = c(5, 3, 2),
                                                                               degree = c(2, 3, 4)),
                                                                   rf = list(num.trees = c(500, 750, 1000),
                                                                             mtry = c(3, 5)),
                                                                   glm = list(alpha = c(1, 0, 0.2, 0.4, 0.6, 0.8)),
                                                                   gbm = list(shrinkage = c(0.005, 0.01, 0.05),
                                                                              n.trees = c(450, 600, 800)),
                                                                   ppr = list(nterms = c(2,3,4),
                                                                              sm.method = c("supsmu")),
                                                                   cubist = list(committees = c(50, 100),
                                                                                 neighbors = c(0, 1, 2, 3)),
                                                                   gp = list(kernel = c("rbfdot", "polydot", "vanilladot"),
                                                                             tol = c(0.01, 0.001)),
                                                                   svm = list(kernel = c("rbfdot","vanilladot"),
                                                                              epsilon = c(0.1, 0.01), C = c(1, 5, 10))),
                                              committee.ratio = .1,
                                              ma.N = 50,
                                              aggregationFUN = "emase-wcommittee"),
                                     Workflow(wf = 'TSE',
                                              learner =  c('MARS', 'PPR', 'RandomForest', 'PPR'),
                                              learner.pars = list(mars = list(nk   = c(5, 3),   degree= c(2, 3, 4)),
                                                                  glm = list(alpha = c(1, 0, 0.2, 0.4, 0.6, 0.8)),
                                                                  rf = list(num.trees = c(500, 750, 1000),
                                                                             mtry = c(3, 5)),
                                                                  ppr = list(nterms = c(2,3,4),
                                                                              sm.method = c("supsmu"))),
                                              aggregationFUN = "emase-scommittee"),
                                     Workflow(wf = 'TSE',
                                              learner = c("baggedtrees"),
                                              learner.pars = list(baggedtrees = list(ntrees = c(1000))))),
                       embedding.dimension = NULL,
                       nReps = 10,
                       modelOutput = TRUE)
  } else if (search.level == "moderate") {
    tsObj <- tsesearch(timeseries = timeseries,
                       workflow =  c(Workflow(wf = 'arimaWorkflow'),
                                     Workflow(wf = 'TSE',
                                              learner =  c('MARS', 'GLM',
                                                           'RandomForest', 'PPR',
                                                           'GBM', 'GP'),
                                              learner.pars <- list(mars = list(nk = c(3, 2),
                                                                               degree = c(3, 4)),
                                                                   rf = list(num.trees = c(750, 1000)),
                                                                   glm = list(alpha = c(1, 0, 0.5)),
                                                                   gbm = list(shrinkage = c(0.005, 0.01, 0.05),
                                                                              n.trees = c(450, 600, 800)),
                                                                   ppr = list(nterms = c(2,3,4),
                                                                              sm.method = c("supsmu")),
                                                                   gp = list(kernel = c("rbfdot", "vanilladot"),
                                                                             tol = c(0.01, 0.001))),
                                              committee.ratio = .1,
                                              ma.N = 50,
                                              aggregationFUN = "emase-wcommittee"),
                                     Workflow(wf = 'TSE',
                                              learner = c("baggedtrees"),
                                              learner.pars = list(baggedtrees = list(ntrees = c(400))))),
                       embedding.dimension = NULL,
                       nReps = 10,
                       modelOutput = TRUE)
  } else {
    tsObj <- tsesearch(timeseries = timeseries,
                       workflow =  c(Workflow(wf = 'arimaWorkflow'),
                                     Workflow(wf = 'TSE',
                                              learner = c("baggedtrees"),
                                              learner.pars = list(baggedtrees = list(ntrees = c(400))))),
                       embedding.dimension = NULL,
                       nReps = 10,
                       modelOutput = TRUE)
  }
  gs.ensembler(tsObj)
}
