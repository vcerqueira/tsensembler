#' out-of-bag function
#'
#' use training data to get oob samples to train the meta-learners
#'
#' @param train train set from the training set
#' @param test test set from the training set
#' @param form formula
#' @param learner base learners
#' @param learner.pars pars
#' @param embedding.dimension k size
#'
#' @return M.ae containing the ae of M and oob.train, which contains
#' the oob samples
#'
#' @export
OOB.fun <- function(train, test, form, learner, learner.pars, embedding.dimension) {
  M <- learnM(form, train, learner, learner.pars, embedding.dimension)
  Y <- get_y(test, form)

  y_hat <- predict(M, test)
  M.ae <- loss_M(y_hat, prop = TRUE, lossFUN = ae)

  list(M.ae = M.ae, oob.train = test)
}

#' oob function for predictions fun
#' it return out of bag predictions
#'
#' @inheritParams OOB.fun
#'
#' @export
OOB.fun.hat <- function(train, test, form, learner, learner.pars, embedding.dimension) {
  M <- learnM(form, train, learner, learner.pars, embedding.dimension)
  Y <- get_y(test, form)

  Y_hat <- prop_hat(predict(M, test))

  cbind.data.frame(target = Y, Y_hat)
}

#' Forward Validation Procedure
#'
#' @param x data to split into \code{nfolds}
#' @param nfolds number of blocks to split data into
#' @param FUN to apply to train/test
#' @param .rbind logical. rbind results from FUN?
#' @param ... further parameters to FUN
#'
#' @export
ForwardValidation <- function(x, nfolds, FUN, .rbind = TRUE, ...) {
  f <- cut(seq_len(NROW(x)), breaks = nfolds, labels = FALSE)

  cv.res <- list()
  seq. <- seq_len(nfolds)
  for (i in seq.[-length(seq.)]) {
    tr.id <- which(f %in% seq_len(i))
    ts.id <- which(f == i + 1L)

    train <- x[tr.id, ]
    test  <- x[ts.id, ]
    cv.res[[i]] <- FUN(train, test, ...)
  }
  if (.rbind) cv.res <- rbind_(cv.res)

  cv.res
}

#' holdout Procedure for oob
#'
#' @param x data
#' @param estimation.ratio ratio of obsrevation to train
#' @param FUN to apply to train/test
#' @param ... further parameters to FUN
#'
#' @export
holdout <- function(x, estimation.ratio, FUN, ...) {
  len <- NROW(x)
  train <- x[ seq_len(estimation.ratio * len), ]
  test <-  x[-seq_len(estimation.ratio * len), ]

  FUN(train, test, ...)
}
