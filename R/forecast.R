#' Combine Base Learners' Predictions
#'
#' This function aggregates the predictions of models according
#' to a given combination type.
#'
#' @param predictions list of predictions from the base models.
#' @param test test set to predict with.
#' @param topmodels list of top learners accross the test set
#' @param modelscores model weights accross the test set
#' @param modelpreweights base learners pre-weights
#' @param aggregationFUN Combination type function
#'
#' @seealso \code{\link{predict}} for details of the predict method for the
#' base models. \code{\link{tseModel-class}} for details how the class of the
#' Time Series Ensemble model.
#'
#' @return final predictions for the test set
#'
#' @export
combinePredictions <- function(predictions,
                               test,
                               topmodels = NULL,
                               modelscores = NULL,
                               modelpreweights = NULL,
                               aggregationFUN) {
  if (!aggregationFUN %in% c("emase-scommittee", "emase-wcommittee",
                           "regret-scommittee", "regret-wcommittee",
                           "emase-all", "regret-all",
                           "static-s", "static-w", "npmase")) {
    stop("Model combination type (\"aggregationFUN\" parameter) is invalid.")
  }
  seq. <- seq_len(NROW(test))

  cpreds <- switch(aggregationFUN,
      "emase-scommittee" = {
        vnapply(seq., function(i) {
          toppreds <- predictions[topmodels[[i]]]
          mean(
            vnapply(toppreds, function(modelpreds) modelpreds[i, "preds"])
          )
        })
      },
      "regret-scommittee" = {
        vnapply(seq., function(i) {
          toppreds <- predictions[topmodels[[i]]]
          mean(
            vnapply(toppreds, function(modelpreds) modelpreds[i, "preds"])
          )
        })
      },
      "emase-all" = {
        vnapply(seq., function(i) {
          testPreds <- vnapply(predictions, function(modelpreds) modelpreds[i, "preds"])
          testWeights <- unlistn(modelscores[i, ])
          sum(testPreds * testWeights)
        })
      },
      "emase-wcommittee" = {
        vnapply(seq., function(i) {
          toppreds <- predictions[topmodels[[i]]]
          testPreds <- vnapply(toppreds, function(modelpreds) modelpreds[i, "preds"])
          testWeights <- proportion(unlistn(modelscores[i, topmodels[[i]]]))
          sum(testPreds * testWeights)
        })
      },
      "npmase" = {
        vnapply(seq., function(i) {
          toppreds <- predictions[topmodels[[i]]]
          testPreds <- vnapply(toppreds, function(modelpreds) modelpreds[i, "preds"])
          testWeights <- proportion(unlistn(modelscores[i, topmodels[[i]]]))
          sum(testPreds * testWeights)
        })
      },
      "regret-all" = {
        vnapply(seq., function(p) {
          testPreds <- vnapply(predictions, function(modelpreds) modelpreds[p, "preds"])
          testWeights <- modelscores[p, ]
          sum(testPreds * testWeights) / sum(testWeights)
        })
      },
      "regret-wcommittee" = {
        vnapply(seq., function(i) {
          testPreds <- vnapply(predictions, function(modelpreds) modelpreds[i, "preds"])
          testWeights <- modelscores[i, topmodels[[i]]]
          sum(testPreds[topmodels[[i]]] * testWeights) / sum(testWeights)
        })
      },
      "static-s" = {
        vnapply(seq., function(i) {
          mean(
            vnapply(predictions, function(modelpreds) modelpreds[i, "preds"])
          )
        })
      },
      "static-w" = {
        vnapply(seq., function(i) {
          testPreds <- vnapply(predictions, function(modelpreds) modelpreds[i, "preds"])
          sum(testPreds * modelpreweights)
        })
      })
  cpreds
}

#' Compute Predictions for a given base model
#'
#' This function computes the predictions of a list of base learners
#' in a test set
#'
#' @param .Model Base learners comprising the ensemble.
#' @param data embedded time series.
#' @param form Formula
#'
#' @return list of predictions for each base learning model.
#'
#' @import glmnet
#' @import ranger
#' @import gbm
#' @import Cubist
#' @import kernlab
#' @export
computePredictions <- function(.Model, data, form) {
  tgt <- get_target(form)
  if (!tgt %in% colnames(data)) stop("Missing target variable on function \"computePredictions\".")
  trues <- data[ ,tgt]

  modelNames <- names(.Model)
  baseModels <- vcapply(modelNames, function(learner) .splitBy_(learner)[1])
  if ("glm" %in% baseModels && is.null(form)) {
    stop("Please provide a valid formula for GLM to make a prediction.")
  }

  preds <- list()
  for (learner in seq_along(.Model)) {
    preds[[learner]] <- switch(baseModels[learner],
      "rf"  = predict(.Model[[learner]], data)$predictions,
      "gbm" = predict(.Model[[learner]], data, n.trees = 100),
      "sae" = deepnet::nn.predict(.Model[[learner]], data),
      "glm" = {
      #if (!"target" %in% colnames(data)) data <- cbind(target = -1, data)
      data <- data[colnames(data) %in% c("target", .Model[[learner]]$beta@Dimnames[[1]])]
      predict.glmnet(.Model[[learner]], model.matrix(form, data), type = "response")
        },
      predict(.Model[[learner]], data)
    )
  }

  preds <- lapply(preds, function(x) {
    aux.df <- data.frame(x, trues)
    colnames(aux.df) <- c("preds", "trues")
    rownames(aux.df) <- NULL
    aux.df
  })
  names(preds) <- names(.Model)

  preds
}
