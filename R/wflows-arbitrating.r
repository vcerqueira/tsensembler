#' arbitrated delegated ensemble for time series forecasting tasks
#'
#' i - learn M
#' ii - get oob samples
#' iii - setup metadata
#' iv - point metalearn
#' v - delegate
#' vi - predict
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
#' @param ... Further parameters to pass to the function
#'
#' @export
boostedADE <- function(form, train, test, learner, learner.pars, ...) {
  K <- get_embedsize(train)
  target <- get_target(form)

  M <- learnM(form, train, learner, learner.pars, K)

  Y <- get_y(test, form)

  Y_hat <- predict(M, test)
  Y_hat.prop <- prop_hat(Y_hat)

  OOB.train <- ForwardValidation(x = train, nfolds = 10, OOB.fun, .rbind = FALSE,
                                 form = form,
                                 learner = learner,
                                 learner.pars = learner.pars,
                                 embedding.dimension = K)

  metadata <- setup.metadata(OOB.train, test, Y_hat)

  oob.train <- rbind_(lapply(OOB.train, function(i) i$oob.train))

  n.oob <- nrow(oob.train)
  n.test <- nrow(test)
  seq.test <- seq_len(n.test)
  timeseq <- seq_len(n.oob + n.test)[-seq_len(n.oob)] - 1 # -1 for gold

  E_hat <- point.metalearn.rf(metadata, timeseq, steps = 10)

  y_hat <- lambda.prediction(test, E_hat, Y_hat.prop, Y, lwindow = 30, lambda = .333)

  res <- list(trues = Y, preds = y_hat)

  res
}

#' arbitrated selecter ensemble
#'
#' selects the most reliable learner
#'
#' @inheritParams boostedADE
#'
#' @export
ASE <- function(form, train, test, learner, learner.pars, ...) {
  K <- get_embedsize(train)
  target <- get_target(form)
  Y <- get_y(test, form)

  M <- learnM(form, train, learner, learner.pars, K)

  Y_hat <- predict(M, test)
  Y_hat.prop <- prop_hat(Y_hat)

  M.ae <- loss_M(Y_hat, prop = TRUE, ae)

  OOB.train <- ForwardValidation(x = train, nfolds = 10, OOB.fun, .rbind = FALSE,
                                 form = form,
                                 learner = learner,
                                 learner.pars = learner.pars,
                                 embedding.dimension = K)

  metadata <- setup.metadata(OOB.train, test, Y_hat, dynamics.FUN = ts.dynamics)

  oob.train <- rbind_(lapply(OOB.train, function(i) i$oob.train))

  n.oob <- nrow(oob.train)
  n.test <- nrow(test)
  seq.test <- seq_len(n.test)
  timeseq <- seq_len(n.oob + n.test)[-seq_len(n.oob)] - 1 # -1 for gold

  E_hat <- point.metalearn.rf(metadata, timeseq, steps = 10)
  W <- t(apply(E_hat, 1, modelWeighting, na.rm = TRUE))
  W <- select_best(W)


  y_hat <- vnapply(seq.test, function(j) sum(Y_hat.prop[j, ] * W[j, ]))

  res <- list(trues = Y, preds = y_hat)

  res
}

#' Forecasting Arbitrated Ensemble
#' no delegation
#'
#' @inheritParams boostedADE
#'
#' @export
AFE <- function(form, train, test, learner, learner.pars, ...) {
  K <- get_embedsize(train)
  target <- get_target(form)
  Y <- get_y(test, form)

  M <- learnM(form, train, learner, learner.pars, K)

  Y_hat <- predict(M, test)
  Y_hat.prop <- prop_hat(Y_hat)

  M.ae <- loss_M(Y_hat, prop = TRUE, ae)

  OOB.train <- ForwardValidation(x = train, nfolds = 10, OOB.fun, .rbind = FALSE,
                                 form = form,
                                 learner = learner,
                                 learner.pars = learner.pars,
                                 embedding.dimension = K)

  metadata <- setup.metadata(OOB.train, test, Y_hat, dynamics.FUN = ts.dynamics)

  oob.train <- rbind_(lapply(OOB.train, function(i) i$oob.train))

  n.oob <- nrow(oob.train)
  n.test <- nrow(test)
  seq.test <- seq_len(n.test)
  timeseq <- seq_len(n.oob + n.test)[-seq_len(n.oob)] - 1 # -1 for gold

  E_hat <- point.metalearn.rf(metadata, timeseq, steps = 10)
  W <- t(apply(E_hat, 1, modelWeighting, na.rm = TRUE))

  y_hat <- vnapply(seq.test, function(j) sum(Y_hat.prop[j, ] * W[j, ]))

  res <- list(trues = Y, preds = y_hat)

  res
}

#' Forecasting Arbitrated Ensemble
#' no pump
#'
#' @inheritParams boostedADE
#'
#' @export
ADE_nopump <- function(form, train, test, learner, learner.pars, ...) {
  K <- get_embedsize(train)
  target <- get_target(form)

  M <- learnM(form, train, learner, learner.pars, K)

  Y_hat <- predict(M, test)
  Y_hat.prop <- prop_hat(Y_hat)
  M.ae <- loss_M(Y_hat, prop = TRUE, ae)
  Y <- get_y(test, form)

  warmup <- 5L

  testcols  <- setdiff(colnames(test), target)
  test.dynamics <- ts.dynamics(test[ ,testcols])
  testset <- cbind.data.frame(test[,testcols], test.dynamics)

  metadata <- lapply(M.ae, function(m) cbind.data.frame(testset, score = m))

  seq. <- seq_len(nrow(test) - 1L)
  seq._ <- seq.[-seq_len(warmup)]

  E_hat <- point.metalearn.rf(metadata, seq._, steps = 10)

  W <- (t(apply(E_hat, 1, modelWeighting, na.rm = TRUE)))
  W <- rbind(matrix(1 / ncol(E_hat),
                    nrow = warmup + 1L,
                    ncol = ncol(E_hat)), W)

  seq.test <- seq_len(nrow(test))
  committee <- meanae.delegation(Y_hat.prop, Y, lambda = 5, committee.ratio = .5)

  y_hat <- vnapply(seq.test, function(j) {
    f.Y_hat <- Y_hat.prop[j , committee[[j]]]
    f.W <- proportion(W[j , committee[[j]]])
    sum(f.Y_hat * f.W)
  })

  res <- list(trues = Y, preds = y_hat)

  res
}

#' Forecasting Arbitrated Ensemble
#' select no pump
#'
#' @inheritParams boostedADE
#'
#' @export
ASE_nopump <- function(form, train, test, learner, learner.pars, ...) {
  K <- get_embedsize(train)
  target <- get_target(form)

  M <- learnM(form, train, learner, learner.pars, K)

  Y_hat <- predict(M, test)
  Y_hat.prop <- prop_hat(Y_hat)
  M.ae <- loss_M(Y_hat, prop = TRUE, ae)
  Y <- get_y(test, form)

  warmup <- 5L

  testcols  <- setdiff(colnames(test), target)
  test.dynamics <- ts.dynamics(test[ ,testcols])
  testset <- cbind.data.frame(test[,testcols], test.dynamics)

  metadata <- lapply(M.ae, function(m) cbind.data.frame(testset, score = m))

  seq. <- seq_len(nrow(test) - 1L)
  seq._ <- seq.[-seq_len(warmup)]

  E_hat <- point.metalearn.rf(metadata, seq._, steps = 10)

  W <- (t(apply(E_hat, 1, modelWeighting, na.rm = TRUE)))
  W <- rbind(matrix(1 / ncol(E_hat),
                    nrow = warmup + 1L,
                    ncol = ncol(E_hat)), W)
  W <- select_best(W)

  seq.test <- seq_len(nrow(test))

  y_hat <- vnapply(seq.test, function(j) sum(Y_hat.prop[j, ] * W[j, ]))

  res <- list(trues = Y, preds = y_hat)

  res
}

#' Forecasting Arbitrated Ensemble
#' simple del
#'
#' @inheritParams boostedADE
#'
#' @export
boostedADE.simpledel <- function(form, train, test, learner, learner.pars, ...) {
  K <- get_embedsize(train)
  target <- get_target(form)

  M <- learnM(form, train, learner, learner.pars, K)

  Y <- get_y(test, form)

  Y_hat <- predict(M, test)
  Y_hat.prop <- prop_hat(Y_hat)

  OOB.train <- ForwardValidation(x = train, nfolds = 10, OOB.fun, .rbind = FALSE,
                                 form = form,
                                 learner = learner,
                                 learner.pars = learner.pars,
                                 embedding.dimension = K)

  metadata <- setup.metadata(OOB.train, test, Y_hat, dynamics.FUN = ts.dynamics)

  oob.train <- rbind_(lapply(OOB.train, function(i) i$oob.train))

  n.oob <- nrow(oob.train)
  n.test <- nrow(test)
  seq.test <- seq_len(n.test)
  timeseq <- seq_len(n.oob + n.test)[-seq_len(n.oob)] - 1 # -1 for gold

  E_hat <- point.metalearn.rf(metadata, timeseq, steps = 10)
  W <- t(apply(E_hat, 1, modelWeighting, na.rm = TRUE))

  committee <- meanae.delegation(Y_hat.prop, Y, lambda = 5, committee.ratio = .5)

  y_hat <- vnapply(seq.test, function(j) {
    f.Y_hat <- Y_hat.prop[j , committee[[j]]]
    f.W <- proportion(W[j , committee[[j]]])
    sum(f.Y_hat * f.W)
  })

  res <- list(trues = Y, preds = y_hat)

  res
}

#' Forecasting Arbitrated Ensemble
#' simple del
#'
#' @inheritParams boostedADE
#'
#' @export
boostedADE.simpledel.nodyns <- function(form, train, test, learner, learner.pars, ...) {
  K <- get_embedsize(train)
  target <- get_target(form)

  M <- learnM(form, train, learner, learner.pars, K)

  Y <- get_y(test, form)

  Y_hat <- predict(M, test)
  Y_hat.prop <- prop_hat(Y_hat)

  OOB.train <- ForwardValidation(x = train, nfolds = 10, OOB.fun, .rbind = FALSE,
                                 form = form,
                                 learner = learner,
                                 learner.pars = learner.pars,
                                 embedding.dimension = K)

  metadata <- setup.metadata(OOB.train, test, Y_hat)

  oob.train <- rbind_(lapply(OOB.train, function(i) i$oob.train))

  n.oob <- nrow(oob.train)
  n.test <- nrow(test)
  seq.test <- seq_len(n.test)
  timeseq <- seq_len(n.oob + n.test)[-seq_len(n.oob)] - 1 # -1 for gold

  E_hat <- point.metalearn.rf(metadata, timeseq, steps = 10)
  W <- t(apply(E_hat, 1, modelWeighting, na.rm = TRUE))

  committee <- meanae.delegation(Y_hat.prop, Y, lambda = 5, committee.ratio = .5)

  y_hat <- vnapply(seq.test, function(j) {
    f.Y_hat <- Y_hat.prop[j , committee[[j]]]
    f.W <- proportion(W[j , committee[[j]]])
    sum(f.Y_hat * f.W)
  })

  res <- list(trues = Y, preds = y_hat)

  res
}


#' Forecasting Arbitrated Ensemble
#' simple del
#'
#' @inheritParams boostedADE
#'
#' @export
boostedADE.simpledel.gp <- function(form, train, test, learner, learner.pars, ...) {
  K <- get_embedsize(train)
  target <- get_target(form)

  M <- learnM(form, train, learner, learner.pars, K)

  Y <- get_y(test, form)

  Y_hat <- predict(M, test)
  Y_hat.prop <- prop_hat(Y_hat)

  OOB.train <- ForwardValidation(x = train, nfolds = 10, OOB.fun, .rbind = FALSE,
                                 form = form,
                                 learner = learner,
                                 learner.pars = learner.pars,
                                 embedding.dimension = K)

  metadata <- setup.metadata(OOB.train, test, Y_hat, dynamics.FUN = ts.dynamics)

  oob.train <- rbind_(lapply(OOB.train, function(i) i$oob.train))

  n.oob <- nrow(oob.train)
  n.test <- nrow(test)
  seq.test <- seq_len(n.test)
  timeseq <- seq_len(n.oob + n.test)[-seq_len(n.oob)] - 1 # -1 for gold

  E_hat <- point.metalearn.gausspr(metadata, timeseq, steps = 10)
  W <- t(apply(E_hat, 1, modelWeighting, na.rm = TRUE))

  committee <- meanae.delegation(Y_hat.prop, Y, lambda = 5, committee.ratio = .5)

  y_hat <- vnapply(seq.test, function(j) {
    f.Y_hat <- Y_hat.prop[j , committee[[j]]]
    f.W <- proportion(W[j , committee[[j]]])
    sum(f.Y_hat * f.W)
  })

  res <- list(trues = Y, preds = y_hat)

  res
}