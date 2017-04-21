#' baggedtrees
#'
#' baggedtrees with different embedding dimensions and summary statistics.
#' This function trains a bagging ensemble of trees using rpart package implementation
#' of decision trees. BaggingDEp_S differs from BaggingE by using different
#' representations of the series. This is accomplished by using
#' different embedding dimensions within the ensemble. With this
#' property the goal is to encourage diversity into the ensemble
#' so to reduce the variance in the combined model.
#' This model is based on the work done by Oliveira and Torgo (2014)
#'
#' @param form formula (e.g. target ~ .).
#' @param data an embedded time series in a data.frame format.
#' @param embedding.dimension size of the embedding vector.
#' @param learner.pars base learner specific set of parameters
#' @param nstats number of stats used to augment the recent time series dynamics information
#' @param ... Further parameters to \code{baggedtrees}
#'
#' @references Oliveira and Torgo (2014) -- Ensembles for Time Series Forecasting.
#' In ACML Proceedings of Asian Conference on Machine Learning. JMLR: Workshop
#' and Conference Proceedings, 2014.
#'
#' @export
#' @import rpart
baggedtrees <- function(form,
                        data,
                        embedding.dimension,
                        nstats,
                        learner.pars, ...) {
  embed.split.at <- c(embedding.dimension, embedding.dimension / 2, embedding.dimension / 4)
  ntrees <- learner.pars[["baggedtrees"]][["ntrees"]]
  n <- NROW(data)
  .i <- 0L
  .n <- floor(ntrees/(2 * length(embed.split.at))); .seqn <- seq_len(.n)
  #.n <- floor(ntrees/(length(embed.split.at))); .seqn <- seq_len(.n)

  .Models <- lapply(embed.split.at, function(K) {
    predictors <- seq_len(K)
    t1 <- lapply(.seqn, function(i) {
      do.call('rpart', c(list(form, data[bootstrap(n), predictors])))
    })
    predictors <- c(seq_len(K),(NCOL(data) - nstats + 1L):NCOL(data))
    t2 <- lapply(.seqn, function(i) {
      do.call('rpart', c(list(form, data[bootstrap(n), predictors])))
    })
    c(t1, t2)
    #t2
  })
  names(.Models) <- paste("decisiontree", seq_along(.Models), sep = "_")

  unlist(.Models, recursive = FALSE)
}

#' SVM
#'
#' SVM is the implementation of learning procedure
#' of an SVM for regression based on the implementation of \code{kernlab} package.
#'
#' @param form formula
#' @param data embedded time series for training
#' @param learner.pars list of learner parameter.
#'
#' @import kernlab
SVM <- function(form, data, learner.pars = NULL) {
  baselearners <- list()
  bl <- 0L
  naming <- NULL

  if (is.null(learner.pars$svm))         learner.pars$svm         <- list()
  if (is.null(learner.pars$svm$scale))   learner.pars$svm$scale   <- FALSE
  if (is.null(learner.pars$svm$type))    learner.pars$svm$type    <- "eps-svr"
  if (is.null(learner.pars$svm$kernel))  learner.pars$svm$kernel  <- "vanilladot"
  if (is.null(learner.pars$svm$epsilon)) learner.pars$svm$epsilon <- 0.1
  if (is.null(learner.pars$svm$C))       learner.pars$svm$C       <- 1

  for (kernel in learner.pars$svm$kernel) {
    for (epsilon in learner.pars$svm$epsilon) {
      for (C in learner.pars$svm$C) {
        bl <- bl + 1L
        baselearners[[bl]] <- do.call("ksvm", c(list(form, data),
                                                scale = learner.pars$svm$scale,
                                                kernel = kernel,
                                                type = learner.pars$svm$type,
                                                epsilon = epsilon,
                                                C = C))
        naming <- c(naming, paste0("svm_", kernel,"_g", epsilon, "c", C))
      }
    }
  }
  names(baselearners) <- naming
  baselearners
}

#' FFNN
#'
#' FFNN is the implementation of learning procedure
#' of an Feed Forward Neural Network for regression, based on the
#' implementation of \code{nnet} package.
#'
#' @param form formula
#' @param data embedded time series
#' @param learner.pars list of learner parameter.
#'
#' @import nnet
FFNN <- function(form, data, learner.pars) {
  baselearners <- list()
  bl <- 0L
  naming <- NULL

  if (is.null(learner.pars$ffnn))        learner.pars$ffnn        <- list()
  if (is.null(learner.pars$ffnn$trace))  learner.pars$ffnn$trace  <- FALSE
  if (is.null(learner.pars$ffnn$linout)) learner.pars$ffnn$linout <- TRUE
  if (is.null(learner.pars$ffnn$size))   learner.pars$ffnn$size   <- 15
  if (is.null(learner.pars$ffnn$decay))  learner.pars$ffnn$decay  <- 0.01
  if (is.null(learner.pars$ffnn$maxit))  learner.pars$ffnn$maxit  <- 1500

  for (maxit in learner.pars$ffnn$maxit) {
    for (size in learner.pars$ffnn$size) {
      for (decay in learner.pars$ffnn$decay) {
        bl <- bl + 1L
        baselearners[[bl]] <- do.call("nnet", c(list(form, data), linout = learner.pars$ffnn$linout, size = size,
                                                maxit = maxit, decay = decay, trace = learner.pars$ffnn$trace, MaxNWts = 1000000))
        naming <- c(naming, paste0("nnet_s", size,"_d", decay, "m", maxit))
      }
    }
  }
  names(baselearners) <- naming
  baselearners
}

#' MARS
#'
#' MARS is the implementation of learning procedure
#' of an Multivariate Adaptive Regression Splines based on the
#' implementation of \code{earth} package.
#'
#' @param form formula
#' @param data embedded time series
#' @param learner.pars list of learner parameter.
#'
#' @import earth
MARS <- function(form, data, learner.pars) {
  baselearners <- list()
  bl <- 0L
  naming <- NULL

  if (is.null(learner.pars$mars))        learner.pars$mars        <- list()
  if (is.null(learner.pars$mars$nk))     learner.pars$mars$nk     <- 15
  if (is.null(learner.pars$mars$degree)) learner.pars$mars$degree <- 3
  if (is.null(learner.pars$mars$thresh)) learner.pars$mars$thresh <- 0.001

  for (nk in learner.pars$mars$nk) {
    for (degree in learner.pars$mars$degree) {
      for (thresh in learner.pars$mars$thresh) {
        bl <- bl + 1L
        baselearners[[bl]] <- do.call("earth", c(list(form, data), nk = nk, degree = degree, thresh = thresh))
        naming <- c(naming, paste0("mars_nk", nk,"_d", degree, "t", thresh))
      }
    }
  }
  names(baselearners) <- naming
  baselearners
}

#' Cubist
#'
#' Cubist is the implementation of learning procedure
#' of an boosted rule-based regression model based on the
#' implementation of \code{Cubist} package.
#'
#' @param form formula
#' @param data embedded time series
#' @param learner.pars list of learner parameter.
#'
#' @import Cubist
Cubist <- function(form, data, learner.pars) {
  baselearners <- list()
  bl <- 0L
  naming <- NULL

  if (is.null(learner.pars$cubist))  learner.pars$cubist     <- list()
  if (is.null(learner.pars$cubist$committees)) learner.pars$cubist$committees <- 30
  if (is.null(learner.pars$cubist$neighbors)) learner.pars$cubist$neighbors <- 0
  form <- stats::as.formula(paste(deparse(form), "-1"))

  predictors <- stats::model.matrix(form, data)
  target <- data[, "target"]
  for (nCommittee in learner.pars$cubist$committees) {
    for (neighbors in learner.pars$cubist$neighbors) {
      bl <- bl + 1L
      baselearners[[bl]] <- do.call("cubist", c(list(predictors, target), committees = nCommittee, neighbors = neighbors))
      naming <- c(naming, paste0("cubist_", nCommittee,"it",neighbors, "nn"))
    }
  }
  names(baselearners) <- naming
  baselearners
}

#' learningRandomForest
#'
#' learningRandomForest is the implementation of learning procedure
#' of an Random Forest model based on the implementation of \code{ranger} package.
#'
#' @param form formula
#' @param data embedded time series
#' @param learner.pars list of learner parameter.
#'
#' @import ranger
RandomForest <- function(form, data, learner.pars) {
  baselearners <- list()
  bl <- 0L
  naming <- NULL

  if (is.null(learner.pars$rf)) learner.pars$rf <- list()
  if (is.null(learner.pars$rf$num.trees)) learner.pars$rf$num.trees <- 100
  if (is.null(learner.pars$rf$mtry))      learner.pars$rf$mtry      <- ceiling(sqrt(NCOL(data)))

  for (num.trees in learner.pars$rf$num.trees) {
    for (mtry in learner.pars$rf$mtry) {
      bl <- bl + 1L
      baselearners[[bl]] <- do.call("ranger", c(list(form, data), num.trees = num.trees, mtry = mtry, write.forest = TRUE))
      naming <- c(naming, paste0("rf_n", num.trees,"_m", mtry))
    }
  }
  names(baselearners) <- naming
  baselearners
}

#' GBM
#'
#' GBM is the implementation of learning procedure
#' of an generalized boosted regression model based on the
#' implementation of \code{gbm} package.
#'
#' @param form formula
#' @param data embedded time series
#' @param learner.pars list of learner parameter.
#'
#' @import gbm
GBM <- function(form, data, learner.pars) {
  baselearners <- list()
  bl <- 0L
  naming <- NULL

  if (is.null(learner.pars$gbm))  learner.pars$gbm <- list()
  if (is.null(learner.pars$gbm$interaction.depth)) learner.pars$gbm$interaction.depth <- 1
  if (is.null(learner.pars$gbm$shrinkage))         learner.pars$gbm$shrinkage <- 0.001
  if (is.null(learner.pars$gbm$n.trees))           learner.pars$gbm$n.trees <- 100

  for (id in learner.pars$gbm$interaction.depth) {
    for (shrinkage in learner.pars$gbm$shrinkage) {
      for (n.trees in learner.pars$gbm$n.trees) {
        bl <- bl + 1L
        baselearners[[bl]] <- do.call("gbm", c(list(form, data), distribution = "gaussian", interaction.depth = id, shrinkage = shrinkage, n.trees = n.trees))
        naming <- c(naming, paste0("gbm_", n.trees,"t", id, "id", shrinkage, "sh"))
      }
    }
  }
  names(baselearners) <- naming
  baselearners
}

#' GLM
#'
#' learningGLM is the implementation of learning procedure
#' of an generalized linear model (e.g. LASSO, Ridge, Elastic-net) based
#' on the implementation of \code{glmnet} package.
#'
#' @param form formula
#' @param data embedded time series
#' @param learner.pars list of learner parameter.
#'
#' @import glmnet
GLM <- function(form, data, learner.pars) {
  data <- soft.completion(data)
  baselearners <- list()
  bl <- 0L
  naming <- NULL

  if (is.null(learner.pars$glm))  learner.pars$glm     <- list()
  if (is.null(learner.pars$glm$alpha)) learner.pars$glm$alpha <- c(0, 1)

  predictors <- model.matrix(target ~., data)
  target <- data[, "target"]

  for (alpha in learner.pars$glm$alpha) {
    bl <- bl + 1L
    min.lambda <- do.call("cv.glmnet", c(list(predictors, target), family = "gaussian", alpha = alpha))$lambda.min
    baselearners[[bl]] <- do.call("glmnet", c(list(predictors, target), alpha = alpha, lambda = min.lambda))
    naming <- c(naming, paste0("glm_", alpha, "alpha"))
  }
  names(baselearners) <- naming
  baselearners
}

#' PPR
#'
#' learningPPR is the implementation of learning procedure
#' of a Projection Pursuit Regression model based on the
#' implementation of \code{stats} package.
#'
#' @param form formula
#' @param data embedded time series
#' @param learner.pars list of learner parameter.
PPR <- function(form, data, learner.pars) {
  baselearners <- list()
  bl <- 0L
  naming <- NULL

  if (is.null(learner.pars$ppr))        learner.pars$ppr <- list()
  if (is.null(learner.pars$ppr$nterms)) learner.pars$ppr$nterms <- 3
  if (is.null(learner.pars$ppr$sm.method)) learner.pars$ppr$sm.method<-"supsmu"
  for (nterm in learner.pars$ppr$nterms) {
    for (smoother in learner.pars$ppr$sm.method) {
      bl <- bl + 1L
      baselearners[[bl]] <- do.call("ppr", c(list(form, data), nterms = nterm, sm.method = smoother))
      naming <- c(naming, paste0("ppr_", nterm,"nterms"))
    }
  }
  names(baselearners) <- naming
  baselearners
}

#' GP
#'
#' learningGP is the implementation of learning procedure
#' of a Gaussian Processes model based on the implementation
#' of \code{kernlab} package.
#'
#' @param form formula
#' @param data embedded time series
#' @param learner.pars list of learner parameter.
#'
#' @import kernlab
GP <- function(form, data, learner.pars) {
  baselearners <- list()
  bl <- 0L
  naming <- NULL

  if (is.null(learner.pars$gp)) learner.pars$gp <- list()
  if (is.null(learner.pars$gp$kernel)) learner.pars$gp$kernel <- "rbfdot"
  if (is.null(learner.pars$gp$tol)) learner.pars$gp$tol <- 0.001
  for (kernel in learner.pars$gp$kernel) {
    for (tolerance in learner.pars$gp$tol) {
      bl <- bl + 1L
      baselearners[[bl]] <- do.call("gausspr", c(list(form, data), type = "regression", kernel = kernel, tol = tolerance))
      naming <- c(naming, paste0("gp_", kernel,"Kernel"))
    }
  }
  names(baselearners) <- naming
  baselearners
}

#' SAE
#'
#' SAE is the implementation of learning procedure
#' of a Stacked Auto Encoder model based on the
#' implementation of \code{deepnet} package.
#' We advise not to use this model unless you are
#' an expert in neural networks. Tunable parameters are
#' Activation function, hidden units, hidden units dropout
#' and visible units dropout.
#'
#' @param form formula
#' @param data embedded time series
#' @param learner.pars list of learner parameter.
#' @param ... Another parameters to pass to deepnet::sae.dnn.train
#'
#' @import deepnet
SAE <- function(form, data, learner.pars) {
  BaseLearners <- list()
  i <- 0L
  naming <- NULL

  if (is.null(learner.pars$sae))  learner.pars$sae     <- list()
  if (is.null(learner.pars$sae$activationfun)) learner.pars$sae$activationfun <- "tanh"
  if (is.null(learner.pars$sae$hidden)) learner.pars$sae$hidden <- c(10)
  if (is.null(learner.pars$sae$hidden_dropout)) learner.pars$sae$hidden_dropout <- 0
  if (is.null(learner.pars$sae$visible_dropout)) learner.pars$sae$visible_dropout <- 0

  predictors <- stats::model.matrix(form, data)
  target <- data[, "target"]
  for (actfun in learner.pars$sae$activationfun) {
    for (hidden.id in seq_along(learner.pars$sae$hidden)) {
      for (hdropout in learner.pars$sae$hidden_dropout) {
        for (vdropout in learner.pars$sae$visible_dropout) {
          bl <- bl + 1L
          BaseLearners[[bl]] <- do.call("sae.dnn.train", c(list(predictors, target), activationfun = actfun,
                                                           hidden = list(learner.pars$sae$hidden[[hidden.id]]),
                                                           hidden_dropout = hdropout, visible_dropout = vdropout, output = "linear"))
          naming <- c(naming, paste0("sae_", actfun,"H",list(learner.pars$sae$hidden[[hidden.id]]), hdropout, vdropout))
        }
      }
    }
  }
  names(BaseLearners) <- vcapply(naming, function(mname) gsub(" ", "", mname))
  BaseLearners
}
