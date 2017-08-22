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
      vcapply(mnames_raw, function(o)
        split_by_(o)[1])

    Y_hat <- list()
    for (bm in seq_along(M)) {
      Y_hat[[bm]] <-
        switch(
          mnames[bm],
          "rf"  = predict(M[[bm]], data)$predictions,
          "gbm" = predict(M[[bm]], data, n.trees = 100),
          "mvr" = {
            predict(M[[bm]])[, , M[[bm]]$best_comp_train]
          },
          "glm" = {
            X_bm <- M[[bm]]$beta@Dimnames[[1]]
            data <- data[colnames(data) %in% c(target_var, X_bm)]
            data_bm <- stats::model.matrix(form, data)

            predict.glmnet(M[[bm]], data_bm, type = "response")
          },
          predict(M[[bm]], data)
        )
    }
    names(Y_hat) <- mnames_raw

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
    seq. <- seq_len(nrow(Y_hat))
    if (!is.null(committee))
      y_hat <-
        vnapply(seq.,
                function(j) {
                  Y_hat_j <- Y_hat[j, committee[[j]]]
                  W_j <- proportion(W[j, committee[[j]]])
                  sum(Y_hat_j * W_j)
                })
    else
      y_hat <-
        vnapply(seq.,
                function(j) {
                  Y_hat_j <- Y_hat[j, ]
                  W_j <- proportion(W[j, ])
                  sum(Y_hat_j * W_j)
                })

    y_hat
  }
