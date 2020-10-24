#' Selecting best model according to weights
#'
#' This function select the best model from a matrix of data x models.
#' For each row (data point), the model with maximum weight is assigned
#' a weight of 1, while the remaining models are assigned a weight of 0.
#'
#' @param model_scores matrix containing the model weights
#' across the observations
#'
#' @family weighting base models
#'
#' @keywords internal
#'
#' @export
select_best <- function(model_scores) {
  rbind_l(
    l1apply(model_scores, function(x) {
      xmax <- which.max(model_scores[x, ])
      model_scores[x, xmax]  <- 1.
      model_scores[x, -xmax] <- 0.
      model_scores[x, ]
    })
  )
}

#' Recent performance of models using EMASE
#'
#' This function computes \strong{EMASE}, Erfc Moving
#' Average Squared Error, to quantify the recent
#' performance of the base models.
#'
#' @param Y_hat A \code{data.frame} containing the predictions of
#' each base model;
#'
#' @param Y know true values from past data to compare the predictions to;
#'
#' @param lambda Window size. Number of periods to average
#' over when computing \strong{MASE};
#'
#' @param omega Ratio of top models in the committee;
#'
#' @param pre_weights The initial weights of the models, computed in
#' the available data during the learning phase;
#'
#' @family weighting base models
#'
#' @return A list containing two objects:
#' \describe{
#' \item{model_scores}{The weights of the models in each time point}
#' \item{top_models}{Models in the committee in each time point}
#' }
#'
#' @export
model_recent_performance <-
  function(Y_hat, Y, lambda, omega, pre_weights) {
    sqr_err <-
      vapply(Y_hat,
             function(p) {
               se(p, Y)
             }, double(NROW(Y_hat)))

    sqr_err <- as.data.frame(sqr_err)

    model_scores <- EMASE(sqr_err, lambda, pre_weights)

    top_models <- get_top_models(model_scores, omega)

    list(model_scores = model_scores,
         top_models = top_models)
  }

#' Weighting Base Models by their Moving Average Squared Error
#'
#' This function computes the weights of the learning models
#' using the Moving Average Squared Error (MASE) function
#' This method provides a simple way to quantify the recent
#' performance of each base learner and adapt the combined
#' model accordingly.
#'
#' @param loss Squared error of the models at each test point;
#'
#' @param lambda Number of periods to average over when computing MASE;
#'
#' @param pre_weights pre-weights of the base models computed in the
#' train set.
#'
#' @return The weights of the models in test time.
#'
#' @family weighting base models
#' @keywords internal
#'
#' @export
EMASE <- function(loss, lambda, pre_weights) {
  MASE <- roll_mean_matrix(loss, lambda)

  scores <-
    apply(MASE,
          1,
          model_weighting, trans = "linear", na.rm = TRUE)

  scores <- data.frame(t(scores))

  scores <- rbind.data.frame(pre_weights, scores[-NROW(scores), ])

  scores
}

#' Extract top learners from their weights
#'
#' This function extracts the top learners at each test point
#' from a score matrix, according to the committee ratio \strong{omega}.
#'
#' @param scores data frame containing the weights;
#' @param omega committee ratio of top base learners
#'
#' @return A list containing the top base models
#' @keywords internal
#'
#' @family weighting base models
#'
#' @export
get_top_models <- function(scores, omega) {
  seq. <- seq_len(nrow(scores))
  threshold <- 1. - omega

  model_scores <-
    transform(scores, beta = apply(scores, 1, stats::quantile, threshold, na.rm = TRUE))

  B <- pmatch("beta", colnames(model_scores))

  top_models <-
    lapply(seq., function(i) {
      x <- model_scores[i, ]
      beta_thresh <- x[["beta"]]
      x["beta"] <- NULL
      unname(which(x >= beta_thresh))
    })

  top_models <-
    lapply(seq_along(top_models),
           function(i) {
             if (length(top_models[[i]]) == 0L)
               top_models[[i]] <- top_models[[i - 1L]]
             top_models[[i]]
           })

  top_models
}


#' Building a committee for an ADE model
#'
#' @param Y_hat A data.frame containing the predictions
#' of base models;
#'
#' @param Y True values of the time interval for
#' which to compute the committee;
#'
#' @param lambda Window size. Number of observations to take into
#' account to build the committee;
#'
#' @param omega Committee ratio -- ratio of models to dynamically weight
#' across the data;
#' @keywords internal
#'
#' @family weighting base models
#'
#' @export
build_committee <-
  function(Y_hat, Y, lambda, omega) {
    Y_loss <- lapply(Y_hat, function(j) ae(j, Y))
    Y_loss <- as.data.frame(Y_loss, row.names = NULL)

    rolled_loss <- roll_mean_matrix(Y_loss, lambda)
    rolled_loss <- rbind.data.frame(rep(1., times = ncol(Y_loss)),
                                    rolled_loss[-NROW(rolled_loss), ])

    beta <- apply(rolled_loss, 1, stats::quantile, probs = omega, na.rm = TRUE)

    C <- l1apply(rolled_loss, function(j) {
      which(rolled_loss[j, ] < beta[j])
    })

    nullC <- which(vapply(C, length, numeric(1)) < 1)
    for (k in nullC)
      C[[k]] <- seq_len(ncol(Y_loss))

    C
  }

#' Model weighting
#'
#' This is an utility function that takes the raw error of
#' models and scales them into a 0-1 range according to one
#' of three strategies:
#'
#' \describe{
#' \item{erfc}{using the complementary Gaussian error function}
#' \item{softmax}{using a softmax function}
#' \item{linear}{A simple normalization using max-min method}
#' }
#'
#' These tranformations culminate into the
#' final weights of the models.
#'
#' @param x A object describing the loss of each base model
#'
#' @param trans Character value describing the transformation type.
#' The available options are \strong{softmax}, \strong{linear} and
#' \strong{erfc}. The softmax and erfc provide a non-linear transformation
#' where the weights decay exponentially as the relative loss of a given model
#' increases (with respect to all available models). The linear transformation
#' is a simple normalization of values using the max-min method.
#'
#' @param ... Further arguments to \code{normalize} and
#' \code{proportion} functions \(na.rm = TRUE\)
#'
#' @return An object describing the weights of models
#'
#' @family weighting base models
#'
#' @export
model_weighting <- function(x, trans = "softmax", ...) {
  if (!trans %in% c("softmax", "linear"))
    stop("Please choose a proper model weighting strategy\n", call. = FALSE)

  if (is.list(x))
    x <- unlist(x, use.names = FALSE)

  if (all(is.na(x))) {
    warning("in model_weighting, all vector is na.")
    return(rep(1/length(x), times = length(x)))
  }

  if (any(is.na(x))) {
    x[is.na(x)] <- max(x, na.rm=TRUE)
  }

  if (trans == "softmax") {
    w <- NA_real_
    while (any(is.na(w))) {
      w <- softmax(-x)
      x <- x / 2
    }
  } else if (trans == "linear") {
    nx <- normalize(-x)
    w <- proportion(nx)
  }

  w[is.na(w)] <- 0.

  w
}

#' Computing the error of base models
#'
#' @param Y_hat predictions of the base models ("@Y_hat" slot)
#' from \code{\link{base_ensemble-class}} object;
#'
#' @param Y true values from the time series;
#'
#' @param lfun loss function to compute. Defaults to \code{ae}, absolute
#' error
#'
#' @keywords internal
#'
#' @export
base_models_loss <-
  function(Y_hat, Y, lfun = se) {

    models_loss <-
        lapply(Y_hat,
               function(o) {
                 lfun(Y, o)
               })

    as.data.frame(models_loss)
  }



#' Get most recent lambda observations
#'
#' @param data time series data as data.frame
#' @param lambda number of observations to keep
#' @keywords internal
#'
#' @export
recent_lambda_observations <-
  function(data, lambda) {
    data[(nrow(data)-lambda+1):nrow(data), ]
  }
