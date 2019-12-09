#' Out-of-bag loss estimations
#'
#' A pipeline for retrieving out-of-bag loss estimations
#'
#' @param train train set from the training set;
#'
#' @param test test set from the training set;
#'
#' @param form formula;
#'
#' @param specs object of class \code{\link{model_specs-class}}. Contains
#' the specifications of the base models.
#'
#' @param lfun loss function for metalearning. Defaults to ae -- absolute error.
#'
#' @param num_cores A numeric value to specify the number of cores used to
#' train base and meta models. num_cores = 1
#' leads to sequential training of models. num_cores > 1
#' splits the training of the base models across num_cores cores.
#'
#' @family out-of-bag functions
#'
#' @keywords internal
#'
#' @return A list containing two objects:
#' \describe{
#' \item{mloss}{loss of base models in \strong{test}}
#' \item{oob}{out-of-bag test samples}
#' \item{Y_hat}{predictions by base models}
#' }
#'
#' @export
intraining_estimations <-
  function(train, test, form, specs, lfun, num_cores) {
    utils::capture.output(M <- build_base_ensemble(form, train, specs, num_cores))
    Y <- get_y(test, form)
    Y_tr <- get_y(train, form)

    Y_hat <- predict(M, test)
    model_loss <-
      base_models_loss(
        Y_hat = Y_hat,
        Y = Y,
        lfun = lfun,
        y_tr = Y_tr
      )

    list(mloss = model_loss,
         oob = test,
         Y_hat = Y_hat)
  }

#' Out-of-bag predictions
#'
#' A pipeline for retrieving out-of-bag predictions
#' from the base models
#'
#' @inheritParams intraining_estimations
#'
#' @family out-of-bag functions
#'
#' @keywords internal
#'
#' @export
intraining_predictions <- function(train, test, form, specs) {
  utils::capture.output(M <- build_base_ensemble(form, train, specs))
  Y <- get_y(test, form)

  Y_hat <- predict(M, test)

  cbind.data.frame(target = Y, Y_hat)
}

#' Prequential Procedure in Blocks
#'
#' @param x data to split into \code{nfolds} blocks;
#'
#' @param nfolds number of blocks to split data into;
#'
#' @param FUN to apply to train/test;
#'
#' @param .rbind logical. If TRUE, the results from
#' FUN are \strong{rbind}ed;
#'
#' @param ... further parameters to FUN
#'
#' @seealso \code{\link{intraining_estimations}}
#' function to use as \strong{FUN} parameter.
#'
#' @keywords internal
#'
#' @export
blocked_prequential <- function(x, nfolds, FUN, .rbind = TRUE, ...) {
  f <- cut(seq_len(NROW(x)), breaks = nfolds, labels = FALSE)

  cv.res <- vector("list", nfolds - 1)
  seq. <- seq_len(nfolds)
  for (i in seq.[-length(seq.)]) {
    tr.id <- which(f %in% seq_len(i))
    ts.id <- which(f == i + 1L)

    train <- x[tr.id, ]
    test  <- x[ts.id, ]
    cv.res[[i]] <- FUN(train, test, ...)
    cv.res[[i]]$test_ids <- ts.id
  }
  if (.rbind) cv.res <- rbind_l(cv.res)

  cv.res
}


#' Holdout
#'
#' @param x data to split into \code{nfolds} blocks;
#' @param beta ratio of observations for training
#' @param FUN function to apply to train/test split
#' @param ... further arguments to FUN
#'
#' @keywords internal
#'
#' @export
holdout <- function(x, beta, FUN, ...) {
  len <- NROW(x)
  train <- x[ seq_len(beta * len), ]
  test <-  x[-seq_len(beta * len), ]

  FUN(train, test, ...)
}
