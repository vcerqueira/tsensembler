#' Dynamic Analysis of Performance of Base Learners
#'
#' This function dynamically weights learning models according to
#' their performance.
#'
#' @param predictions A list of data.frames comprising the predictions
#' of the base models. Each \code{data.frame} contains
#' information on the performance of a model: One column with the true values
#' and a second column with the respective predictions
#' @param ma.N Number of periods to average over when computing MASE.
#' @param committee.ratio Proportion of models (top performing ones) used in
#' each prediction time.
#' @param preweights The initial weights of the models, computed in
#' the available data.
#' @param aggregationFUN Manner of aggregation of learners. Basic form is to
#' avegare the prediction of all learners.
#' @param Mstar todo
#'
#' @return The weights of each learning model (fScores). If applicable, it also returns
#' the top performing models (TOPmodels)
#'
#' @export
LearnerAnalysis <- function(predictions,
                            ma.N,
                            committee.ratio = .1,
                            preweights,
                            aggregationFUN,
                            Mstar = NULL) {
  if (committee.ratio <= 0 || committee.ratio > 1) {
    stop("Please choose a valid top committee ratio (double between 0 and 1).")
  }

  SE <- as.data.frame(vapply(predictions, function(p) {
    se(p[["trues"]], p[["preds"]])
  }, double(NROW(predictions[[1]]))))

  if (aggregationFUN %in% c("regret-scommittee", "regret-all", "regret-wcommittee")) {
    fScores <- Regret(SE, Mstar)
  } else if (aggregationFUN %in% c("emase-scommittee", "emase-wcommittee", "emase-all")) {
    fScores <- EMASE(SE, ma.N, preweights)
  } else if (aggregationFUN %in% c("npmase")) {
    fScores <- NP_MASE(SE, ma.N, preweights)
  }

  Toplearners <- getToplearners(fScores, committee.ratio)

  list(Toplearners = Toplearners, fScores = fScores)
}

#' Weighting Base Models by their Moving Average Squared Error
#'
#' This function computes the weights of the learning models
#' using the Moving Average Squared Error \(MASE\) function
#' This method provides a simple way to quantify the recent
#' performance of each base learner and adapt the combined
#' model accordingly.
#'
#' @param SE Squared Error of the models at each test point
#' @param ma.N Number of periods to average over when computing MASE
#' @param preweights pre-weights of the base models computed in the
#' train set.
#' @param smartass if TRUE, lastpoint is used to compute EMASE, 
#' breaking the golden rule. Defaults to FALSE. Useful for hindsight analysis.
#'
#' @return The weights of the models in test time.
#'
#' @export
EMASE <- function(SE, ma.N, preweights, smartass = FALSE) {
  MASE <- rollmeanmatrix(SE, ma.N)

  fScores <- data.frame(t(apply(MASE, 1, modelWeighting)))
  if (!smartass) {
    fScores <- rbind.data.frame(preweights, fScores[-NROW(fScores), ])
  }

  fScores
}

#' Weighting Base Models by their Moving Average Squared Error
#'
#' This function computes the weights of the learning models
#' using the Moving Average Squared Error \(MASE\) function
#' This method provides a simple way to quantify the recent
#' performance of each base learner and adapt the combined
#' model accordingly.
#'
#' @param SE Squared Error of the models at each test point
#' @param ma.N Number of periods to average over when computing MASE
#' @param preweights pre-weights of the base models computed in the
#' train set.
#' @param smartass if TRUE, lastpoint is used to compute EMASE, 
#' breaking the golden rule. Defaults to FALSE. Useful for hindsight analysis.
#'
#' @return The weights of the models in test time.
#'
#' @export
NP_MASE <- function(SE, ma.N, preweights, smartass = FALSE) {
  MASE <- rollmeanmatrix(SE, ma.N)

  fScores <- data.frame(t(apply(MASE, 1, function(j) {
    proportion(normalizeMaxMin(-j))
  })))
  if (!smartass) {
    fScores <- rbind.data.frame(preweights, fScores[-NROW(fScores), ])
  }
  fScores
}

#' Weighting base Models by Regret
#'
#' @param SE Squared Error of the models at each test point
#' @param Mstar Maximum expected loss
Regret <- function(SE, Mstar) {
  SE[ ] <- lapply(SE / Mstar, cumsum)
  SE   <- rbind.data.frame(rep(1., times = ncol(SE)), SE[-nrow(SE), ])
  N <- ncol(SE)
  seq. <- seq_len(NROW(SE))

  fScores <- as.data.frame(t(vapply(seq., function(i) {
    iloss <- unlistn(SE[i, ])
    expLoss(iloss, N, i, .proportion = TRUE)
  }, double(NCOL(SE)))))
  colnames(fScores) <- colnames(SE)

  fScores
}

#' Extract Top Learners from the Weights of Base Models
#'
#' This function extracts the top learners at each test point
#' from a score matrix, according to the committee ratio.
#'
#' @param scores data frame containing the scores
#' @param committee.ratio ratio of base learners to keep
#'
#' @return Top learners in test time
#' @export
getToplearners <- function(scores, committee.ratio) {
  seq. <- seq_len(nrow(scores))
  threshold <- 1. - committee.ratio
  model_scores <- transform(scores, beta = apply(scores, 1, quantile, threshold))
  B <- pmatch("beta", colnames(model_scores))

  Toplearners <- lapply(seq., function(i) {
    x <- model_scores[i, ]
    beta_thresh <- x[["beta"]]
    x["beta"] <- NULL
    unname(which(x >= beta_thresh))
  })

  Toplearners <- lapply(seq_along(Toplearners), function(i) {
    if (length(Toplearners[[i]]) == 0L) Toplearners[[i]] <- Toplearners[[i - 1L]]
    Toplearners[[i]]
  })

  Toplearners
}

#' Learning the Base Learners that comprise the Ensemble
#'
#' Learntse function applies the learning algorithms
#' inside the \code{TSE} workflow. It receives as input a training dataset
#' and the meta-data about the learning algorithms. The output is a series of
#' predictive models, the relative performance of these models
#' in the training data and the maximum expected loss incurred by
#' the learning models in the data.
#'
#' @param train Training set
#' @param form Formula
#' @param learner Character describing the learning models
#' @param learner.pars List of parameters regarding the learning models.
#' @param Ksplits Column-wise (predictors) splits to be made
#' @param Rsplits Row-wise (training window) splits to be made
#' @param nstats Number of descriptive statistics used.
#' @param embedding.dimension todo
#' @param verbose Logical. If TRUE, some information about the current status of the
#' training procedure will be printed into the console.
#'
#' @seealso \code{\link{tseModel-class}} on detail for building an ensemble
#' with the results from \emph{Learntse}. \code{\link{predict}} for details on
#' the predict method for \emph{tseModel} class objects.
#'
#' @return A series of predictive models (\code{learningModels}), the weights of
#' the models computed in the training data (\code{preweights}) and the
#' maximum expected loss of all models in the data (\code{Mstar})
#' @export
Learntse <- function(train,
                     form,
                     learner,
                     learner.pars,
                     Ksplits = NULL,
                     Rsplits = NULL,
                     nstats,
                     embedding.dimension,
                     verbose) {
  if (learner %in% "baggedtrees" && length(learner) == 1L) {
      .Models <- do.call(learner, list(form,
                                   train,
                                   embedding.dimension = embedding.dimension,
                                   nstats = nstats,
                                   learner.pars = learner.pars,
                                   cp = 0,
                                   minsplit = 6))
    return(list(learningModels = .Models, preweights = NULL, Mstar = NULL))
  }
  if (is.null(Rsplits)) Rsplits <- nrow(train)
  if (is.null(Ksplits)) Ksplits <- embedding.dimension
  .Models <- list()
  .tPerf  <- list()
  .i <- 0L
  total.models <- length(learner) * length(Ksplits) * length(Rsplits)
  row.size <- NROW(train)
  col.size <- NCOL(train)

  summary.stats.id <- ifelse(nstats > 0, c(col.size, col.size - 1L), integer(0L))
  for (K in Ksplits) {
    for (R in Rsplits) {
      for (model in learner) {
        .i <- .i + 1L
        if (verbose) {
          cat("Training model no. ", .i,
              " of a total of ", total.models, ".\n")
        }
        .rows <- (row.size - R):row.size
        .cols <- na.omit(c(seq_len(K), summary.stats.id))
        .Models[[.i]] <- do.call(model, c(list(form, train[.rows, .cols], learner.pars)))
        .tPerf[[.i]] <- computePredictions(.Models[[.i]], train[.rows, .cols], form)
      }
    }
  }
  .Models <- unlist(.Models, recursive = FALSE)
  .tPerf  <- unlist(.tPerf, recursive = FALSE)

  Mstar <- Reduce(max,
                  lapply(.tPerf, function(x) {
                    ae(true = x[,'trues'], pred = x[,'preds'])
                    })
                  )

  prew <- modelWeighting(vnapply(.tPerf,
                                function(r) mse(r[["trues"]], r[["preds"]])))

  list(learningModels = .Models, preweights = prew, Mstar = Mstar)
}

#' Selecting best model according to weights
#'
#' @param model_scores Matrix containing the model weights across the data
#'
#' @export
select_best <- function(model_scores) {
  rbind_(
    l1apply(model_scores, function(x) {
      xmax <- which.max(model_scores[x, ])
      model_scores[x, xmax]  <- 1.
      model_scores[x, -xmax] <- 0.  
      model_scores[x, ]
    })
  )
}


#' Wrapper for learning the base learning models M
#' It encapsulates the functions \code{Learntse} and \code{tseModel}
#'
#' @param form formula
#' @param data training data
#' @param learner learner info
#' @param learner.pars parameter setting for the models
#' @param embedding.dimension k size
#'
#' @export
learnM <- function(form, data, learner, learner.pars, embedding.dimension) {
  Learningensemble <- Learntse(train = data,
                               form = form,
                               learner = learner,
                               learner.pars = learner.pars,
                               nstats = 0,
                               embedding.dimension = embedding.dimension,
                               verbose = FALSE)

  tseModel(baseModels = Learningensemble$learningModels,
           preWeights = Learningensemble$preweights,
           form = form,
           ma.N = NULL,
           embedding.dimension = embedding.dimension,
           committee.ratio = NULL,
           aggregationFUN = "static-s",
           Mstar = Learningensemble$Mstar)
}

#' from Y_hat and Y calculate mean ae and committee
#'
#' @param Y_hat predictions
#' @param Y true values
#' @param lambda number of observations
#' @param committee.ratio committee ratio
#'
#' @export
meanae.delegation <- function(Y_hat, Y, lambda, committee.ratio = .5) {
  ae.M <- as.data.frame(lapply(Y_hat, function(j) ae(j, Y)))
  rownames(ae.M) <- NULL
  ae.M.smooth <- rollmeanmatrix(ae.M, lambda)
  ae.M.smooth <- rbind.data.frame(rep(1., times = ncol(ae.M)), 
                                  ae.M.smooth[-NROW(ae.M.smooth), ])
  
  beta <- apply(ae.M.smooth, 1, quantile, probs = committee.ratio)
  
  C <- l1apply(ae.M.smooth, function(j) {
    which(ae.M.smooth[j, ] < beta[j])
  })
  
  nullC <- which(vnapply(C, length) < 1)
  for (k in nullC) C[[k]] <- seq_len(ncol(ae.M))
  
  C
}

#' Model Weighting with Erfc
#'
#' This is an utility function that takes the raw error of models and scales
#' them into a 0-1 range. 
#' 
#' One can use a simple linear transformation (\code{linear}) or more
#' non-linear ones such as the \code{erfc} function or the \code{softmax}. 
#' These tranformations culminate into the final weights of the models.
#' 
#' @param x Loss metric value
#' 
#' @param trans Character value describing the transformation type. 
#' The available options are \strong{softmax}, \strong{linear} and
#' \strong{erfc}. The softmax and erfc provide a non-linear transformation
#' where the weights decay exponentially as the relative loss of a given model 
#' increases (with respect to all available models). The linear transformation
#' is a simple normalization of values.
#' 
#' @param ... Further arguments to \code{normalizeMaxMin} and
#' \code{proportion} functions (e.g. na.rm = TRUE)
#'
#' @return Weights of models
#' 
#' @export
model_weighting <- function(x, trans = "softmax", ...) {
  if (!trans %in% c("softmax", "linear", "erfc")) 
    stop("Please choose a proper model weighting strategy\n", call. = FALSE)
  
  if (is.list(x)) 
    x <- unlistn(x)
  
  if (trans == "softmax") {
    w <- softmax(-x)
  } else if (trans == "linear") {
    nx <- normalizeMaxMin(-x, ...)
    w <- proportion(nx)
  } else {
    nx <- normalizeMaxMin(x, ...)
    e_nx <- erfc(nx)
    w <- proportion(e_nx, ...)
  }
  
  w[is.na(w)] <- 0.
  
  w
}
