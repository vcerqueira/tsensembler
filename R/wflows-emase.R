#' Time Series Ensemble
#'
#' This function comprises the pipeline used to construct the
#' time series ensemble.
#'
#' @param form Formula
#' @param train embedded time series used for training the base learners
#' @param test embedded time series used for testing
#' @param learner Character vector describing the base algorithms to be trained.
#' Current available implemented models are:
#' \describe{
#'    \item{MARS}{Multivariate Adaptive Regression Splines from \strong{earth} package}
#'    \item{PPR}{Projection Pursuit Regression from \strong{stats} package}
#'    \item{baggedtrees}{For a bagging ensemble for time series forecasting tasks
#'     (Oliveira and Torgo, 2014)}
#'    \item{SVM}{Support Vector Machines from \strong{kernlab} package}
#'    \item{GP}{Gaussian Processes from \strong{kernlab} package}
#'    \item{FFNN}{Feed Forward Neural Networks from \strong{nnet} package}
#'    \item{SVM}{Support Vector Machines from \strong{kernlab} package}
#'    \item{Cubist}{Rule-based Regression from \strong{Cubist} package}
#'    \item{RandomForest}{Random Forests from \strong{ranger} package}
#'    \item{GBM}{Generalized Boosted Regression from \strong{gbm} package}
#'    \item{GLM}{Generalized Linear Models (e.g. Ridge regression, LASSO, Elastic-Net)
#'    from \strong{glmnet} package}
#'    \item{SAE}{Stacked Autoencoder from \strong{deepnet} package}
#' }
#' @param learner.pars Named list describing the parameter of the \code{learner}. Below are
#' described some examples.
#' @param varying.embed Logical: Should different embedding dimensions be used. Defaults to TRUE.
#' @param varying.trainwindow Logical: Should varying training windows be used. Defaults to TRUE.
#' @param ma.N Integer: Number of periods to average over when computing MASE.
#' @param committee.ratio A numeric value between 0 and 1 representing the ratio
#' of base learners that comprise the committee at each prediction time.
#' If \code{committee.ratio} equals, say, 0.2, at time \emph{t}, the ensemble will the
#' 20\% best base learners up to time \emph{t - 1}.
#' @param aggregationFUN The function name used to combine the base learners. See
#' @param modelOutput Logical. If TRUE, the methods applied to \code{tsensembler-class} return
#' a predictive model ensemble. Otherwise, it returns performance statistics based on
#' \code{\link[performanceEstimation]{performanceEstimation}}.
#' @param ... Further parameters to pass to \code{TSE}
#' @param verbose If TRUE, the status of the training of base learners is printed to the console.
#'
#' @seealso \code{\link{ensembler}} for details on the method that encapsulates this
#' workflow.
#'
#' @examples
#' \dontrun{
#' workflow <- Workflow(wf = 'TSE',
#'                     learner = c('MARS', 'PPR'),
#'                     learner.pars = list(mars = list(nk = c(5, 2),
#'                                                     degree= c(3, 4)),
#'                                         ppr = list(nterms = c(2,3,4))),
#'                     varying.embed = TRUE,
#'                     varying.trainwindow = FALSE,
#'                     committee.ratio = .1,
#'                     aggregationFUN = "regret",
#'                     verbose = FALSE)
#' }
#'
#' @return a list containing the true values and predicted values for test set
#'
#' @export
TSE <- function(form,
                train,
                test,
                learner = NULL,
                learner.pars = NULL,
                varying.embed = TRUE,
                varying.trainwindow = TRUE,
                ma.N = NULL,
                committee.ratio = NULL,
                aggregationFUN,
                verbose = TRUE,
                modelOutput = FALSE, ...) {
  if (!is.null(aggregationFUN) && !aggregationFUN %in% c("emase-scommittee", "emase-wcommittee",
                                  "regret-scommittee", "regret-wcommittee",
                                  "emase-all", "regret-all",
                                  "static-s", "static-w", "npmase")) {
    stop("Model combination type (\"aggregationFUN\" parameter) is invalid.")
  }
  tgt <- get_target(form)

  if (is.null(learner.pars)) learner.pars <- list()
  if (!methods::is(learner.pars, "list")) learner.pars <- as.list(learner.pars)

  embed.cols <- which(get_embedcols(train))

  embedding.dimension <- length(embed.cols) + 1
  row.size   <- NROW(train)

  traintarget <- train[ ,tgt]
  testtarget  <- test[ ,tgt]

  train <- cbind.data.frame(target = traintarget, embedStats(train[ ,embed.cols]))
  test  <- cbind.data.frame(target = testtarget,  embedStats(test[ ,embed.cols]))

  nstats <- ncol(train) - embedding.dimension

  trues <- model.response(model.frame(form, test, na.action = NULL))

  if (varying.embed) {
    Ksplits <- c(embedding.dimension, embedding.dimension/2, embedding.dimension/3)
  } else {
    Ksplits <- embedding.dimension
  }
  if (varying.trainwindow) {
    Rsplits <- floor(c(row.size, row.size/2, row.size/4)) - 1
  } else {
    Rsplits <- row.size
  }

  Learningensemble <- Learntse(train,
                               form,
                               learner,
                               learner.pars,
                               Ksplits,
                               Rsplits,
                               nstats,
                               embedding.dimension,
                               verbose)

  .tse <- tseModel(baseModels = Learningensemble$learningModels,
                   preWeights = Learningensemble$preweights,
                   form = form,
                   ma.N = ma.N,
                   embedding.dimension = embedding.dimension,
                   committee.ratio = committee.ratio,
                   aggregationFUN = aggregationFUN,
                   Mstar = Learningensemble$Mstar)

  predictions <- predict(.tse, test)

  res <- list(trues = trues, preds = predictions@y_hat)
  modelOutput <- FALSE
  if (modelOutput) res <- c(res, tse = .tse)

  res
}


