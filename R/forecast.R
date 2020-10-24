#' Compute the predictions of base models
#'
#' This function is used to predict new observations
#' using the predictive models comprising an ensemble.
#' It calls on the respective method based on the type of
#' model, and returns the predictions as a list.
#'
#' @param M list of base models;
#' @param form formula;
#' @param data new data to predict;
#'
#' @import glmnet
#' @import ranger
#' @import gbm
#' @import Cubist
#' @import kernlab
#' @import pls
#'
#' @keywords internal
#'
#' @export
compute_predictions <-
  function(M, form, data) {

    target_var <- get_target(form)

    mnames_raw <- names(M)
    mnames <-
      vapply(mnames_raw,
             function(o) {
               split_by_(o)[1]
             }, character(1), USE.NAMES = FALSE)

    Y_hat <- list()
    time_cost <- c()
    for (bm in seq_along(M)) {
      t0 <- Sys.time()
      Y_hat[[bm]] <-
        switch(
          mnames[bm],
          "rf"  = predict(M[[bm]], data)$predictions,
          "gbm" = predict(M[[bm]], data, n.trees = 100),
          "mvr" = {
            predict_pls_pcr(M[[bm]], data)
          },
          "nnet" = {
            newX <- as.matrix(stats::model.matrix(form, data))
            monmlp.predict(newX, M[[bm]])[,1]
          },
          "xgb" = {
            xgb_predict(M[[bm]], data)
          },
          "glm" = {
            X_bm <- M[[bm]]$beta@Dimnames[[1]]
            data <- data[colnames(data) %in% c(target_var, X_bm)]
            data_bm <- stats::model.matrix(form, data)

            predict.glmnet(M[[bm]], data_bm, type = "response")
          },

          "cub" = {
            data_bm <- stats::model.matrix(form, data)
            predict(M[[bm]], data_bm, neighbors=M[[bm]]$neighbors)
          },
          predict(M[[bm]], data)
        )

      t1 <- Sys.time()
      dt <- difftime(t1,t0, units="secs")
      dt <- round(as.vector(dt),5)

      time_cost <- c(time_cost,dt)
    }
    names(Y_hat) <- mnames_raw
    names(time_cost) <- mnames_raw

    attr(Y_hat, "Times") <- time_cost

    Y_hat
  }

#' Combining the predictions of several models
#'
#' This function simply applies a weighted average,
#' where the predictions of the base models \strong{Y_hat}
#' are weighted according to their weights \strong{W}. If
#' a \strong{committee} is specified, only models in the committee are
#' weighted.
#'
#' @param Y_hat a data.frame with the predictions of the
#' base models;
#'
#' @param W a matrix or data.frame with the weights of the
#' base models;
#'
#' @param committee A list containing the ids of the models in the
#' committee.
#'
#' @keywords internal
#'
#' @export
combine_predictions <-
  function(Y_hat, W, committee = NULL) {
    if (nrow(Y_hat) != nrow(W)) {
      stop("in prediction comb, nrow YH differs from W")
    }

    if (!is.null(committee) && length(committee) != nrow(W)) {
      stop("in pred comb, lengh C diffs nrow W")
    }


    seq. <- seq_len(nrow(Y_hat))
    if (!is.null(committee))
      y_hat <-
        vapply(seq.,
                function(j) {
                  Y_hat_j <- Y_hat[j, committee[[j]]]
                  W_j <- proportion(W[j, committee[[j]]])
                  sum(Y_hat_j * W_j)
                }, numeric(1), USE.NAMES = FALSE)
    else
      y_hat <-
        vapply(seq.,
                function(j) {
                  Y_hat_j <- Y_hat[j, ]
                  W_j <- proportion(W[j, ])
                  sum(Y_hat_j * W_j)
                }, numeric(1), USE.NAMES = FALSE)

    y_hat
  }
