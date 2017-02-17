#' setup metadata using oob samples and ae of models
#'
#' using forward validation we get the oob train data and
#' oob ae of models to augment our meta set
#' here we combine these to the test info to setup the meta data
#' to train the overlineM
#'
#' @param OOB.train oob train and ae info returned by \code{ForwardValidation} fun,
#' combined with \code{OOB.fun}.
#' @param test test set
#' @param Y_hat M prediction for the test set
#' @param dynamics.FUN function that calculates dynamics
#'
#' @export
#'
#' @return metadata ready for metalearning the overlineM
setup.metadata <- function(OOB.train, test, Y_hat, dynamics.FUN = NULL) {
  Y_hat.ae <- loss_M(Y_hat, prop = TRUE, ae)

  oob.train <- rbind_(lapply(OOB.train, function(i) i$oob.train))
  oob.Y_hat.ae <- rbind_(lapply(OOB.train, function(i) i$M.ae))

  ext.test <- rbind.data.frame(oob.train, test)
  
  if (!is.null(dynamics.FUN)) {
    cols <- get_embedcols(ext.test)

    dynamics <- dynamics.FUN(ext.test[ ,cols])
    ext.test <- cbind.data.frame(ext.test, dynamics)
  }

  ext.L_hat <- as.list(rbind.data.frame(oob.Y_hat.ae, Y_hat.ae))

  ext.test$target <- NULL

  metadata <- lapply(ext.L_hat, function(l_hat) cbind.data.frame(ext.test, score = l_hat))

  metadata
}


#' lambda prediction of arbitrated forecasters
#'
#' the idea is to track the performance of the base models
#' those who have an accuracy - which is defined as REF - below a
#' lambda threshold are left out of the prediction.
#'
#' @param test test set
#' @param E_hat oM predictions of M ae
#' @param Y_hat predictions made by the base models
#' @param Y true values
#' @param lwindow recent window size to compute \strong{lambda threshold}.
#' @param lambda accuracy threshold
#'
#' @export
lambda.prediction <- function(test, E_hat, Y_hat, Y, lwindow, lambda) {
  seq.test <- seq_len(nrow(test))

  W <- t(apply(E_hat, 1, modelWeighting, na.rm = TRUE))
  delegated.Y_hat <- lambda.delegation(Y_hat, Y, test, lwindow, lambda)

  vnapply(seq.test, function(j) {
    f.Y_hat <- Y_hat[j , delegated.Y_hat[[j]]]
    f.W <- proportion(W[j , delegated.Y_hat[[j]]])
    sum(f.Y_hat * f.W)
  })
}

#' lambda delegation function
#'
#' @param Y_hat predictions made by the base models
#' @param Y true values
#' @param test test set
#' @param lwindow recent window size to compute \strong{lambda threshold}.
#' @param lambda accuracy threshold
#'
#' @export
lambda.delegation <- function(Y_hat, Y, test, lwindow, lambda) {
  seq.test <- seq_len(nrow(test))

  yhat_acc <- model_accuracy(Y_hat, Y)

  rollacc <- rollmeanmatrix(yhat_acc, lwindow)
  rollacc <- rbind.data.frame(rep(1, times = ncol(rollacc)), rollacc[-nrow(rollacc), ])

  lapply(seq.test, function(y) {
    apt_yy <- which(rollacc[y, ] > lambda)
    if (length(apt_yy) > 0)
      res <- apt_yy
    else
      res <- seq_along(rollacc)

    res
  })
}

