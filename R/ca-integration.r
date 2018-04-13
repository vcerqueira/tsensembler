#' Get predict data for generalising
#'
#' @param object a \code{\link{constructive_aggregation-class}} object
#' @param newdata new data for prediction
#'
#' @keywords internal
#'
#' @export
hat_info <-
  function(object, newdata) {
    Y_hat <- predict(object@base_ensemble, newdata)
    Y <- get_y(newdata, object@form)

    Y_hat_recent <- predict(object@base_ensemble, object@recent_series)
    Y_recent <- get_y(object@recent_series, object@form)

    Y_hat_ext <- rbind.data.frame(Y_hat_recent, Y_hat)
    Y_ext <- c(Y_recent, Y)

    E_ext <- base_models_loss(Y_hat_ext, Y_ext, ae)

    W_E <-
      window_loss(loss = E_ext,
                  lambda = object@lambda,
                  w_0 = object@base_ensemble@pre_weights,
                  trans = "linear")

    C_hat <-
      merging_in_experts(C = object@committee_set$C,
                         Y_hat = Y_hat_ext,
                         W = W_E,
                         FUN = object@aggregate_subsets)

    colnames(C_hat) <- paste0("C_", seq_len(ncol(C_hat)))


    E_C <- base_models_loss(C_hat, Y_ext, ae)
    E_0 <- rep(1 / ncol(E_C), times = ncol(E_C))
    W_C <-
      window_loss(loss = E_C,
                  lambda = object@lambda,
                  w_0 = E_0,
                  trans = "linear")

    C_hat <- utils::tail(C_hat, nrow(C_hat) - object@lambda)
    W_C   <- utils::tail(W_C, nrow(W_C) - object@lambda)
    W_E   <- utils::tail(W_E, nrow(W_E) - object@lambda)

    list(C_hat = C_hat, W_C = W_C, W_E = W_E, Y_hat = Y_hat, Y=Y)
  }

#' Merge models in each committee
#'
#' @param C list of ids representing
#' the models in each subset
#' @param Y_hat data frame with the prediction of each model
#' across data
#' @param W weights of each model according to recent performance
#' @param FUN function used to combine the models in each
#' subset. Currently implemented options are: window_loss,
#' simple, and blast.
#'
#' @keywords internal
#'
#' @export
merging_in_experts <-
  function(C, Y_hat, W, FUN) {
    C <-
      switch(FUN,
             "window_loss" = {
               lapply(C,
                      function(x) {
                        C_ext <-
                          replicate(nrow(W),
                                    x,
                                    simplify = FALSE)

                        combine_predictions(Y_hat, W, C_ext)
                      })
             },
             "blast" = {
               lapply(C,
                      function(x) {
                        W_c <- W[,x]
                        C_hat <- Y_hat[,x]

                        W_c <- select_best(W_c)

                        combine_predictions(C_hat, W_c, NULL)
                      })
             },
             "simple" = {
               lapply(C,
                      function(x) {
                        if (length(x) < 2) {
                          Y_hat[,x]
                        } else {
                          rowMeans(Y_hat[,x])
                        }
                      })
             })

    C <- as.data.frame(C)
    colnames(C) <- paste0("C_", seq_len(ncol(C)))

    C
  }

#' Merge across sub-ensembles
#'
#' @param predict_info data from ...
#' @param FUN function used to combine the models in each
#' subset.
#' @param object constructive aggregation class object
#' @param newdata new data to make predictions
#'
#' @keywords internal
#'
#' @export
combine_committees <-
  function(predict_info, FUN, object, newdata) {
    available_funs <-
      c("window_loss", "simple","blast",
        "arbitrage",
        "loss_train","oracle","mlpol","ewa",
        "fs","ogd","ridge")

    if (!FUN %in% available_funs) {
      stop("combiner not found.")
    }

    y_hat <-
      switch(FUN,
             "window_loss" = {
               combine_predictions(Y_hat = predict_info$C_hat,
                                   W = predict_info$W_C,
                                   committee = NULL)
             },
             "blast" = {
               W_C <- select_best(predict_info$W_C)

               combine_predictions(Y_hat = predict_info$C_hat,
                                   W = W_C,
                                   committee = NULL)
             },
             "simple" = {
               rowMeans(predict_info$C_hat)
             },
             "arbitrage" = {
               CA.ADE_hat(object, predict_info, newdata)
             },
             "mlpol" = {
               CA.MLpol_hat(object, predict_info)
             },
             "ridge" = {
               CA.Ridge_hat(object, predict_info)
             },
             "ogd" = {
               CA.OGD_hat(object, predict_info)
             },
             "ewa" = {
               CA.EWA_hat(object, predict_info)
             },
             "fs" = {
               CA.FixedShare_hat(object, predict_info)
             })

    y_hat
  }

window_loss <-
  function(loss, lambda, w_0, trans = "linear") {
    if (ncol(loss) == 1) {
      coln <- colnames(loss)
      aux_df<-data.frame(aux=rep(1, times=nrow(loss)))
      colnames(aux_df)<- coln
      return(aux_df)
    }

    MASE <- roll_mean_matrix(loss, lambda)

    W <-
      apply(MASE,
            1,
            model_weighting, trans = trans, na.rm = TRUE)

    W <- data.frame(t(W))

    W <- rbind.data.frame(w_0, W[-nrow(W),])

    W
  }

OOB_Chat <-
  function(object) {

    OOB_Y_hat <- object@out_of_bag$Y_hat
    OOB_E <- object@out_of_bag$mloss
    OOB_W <-
      window_loss(loss = OOB_E,
                  lambda = object@lambda,
                  w_0 = object@base_ensemble@pre_weights,
                  trans = "linear")

    OOB_C_hat <-
      merging_in_experts(C = object@committee_set,
                         Y_hat = OOB_Y_hat,
                         W = OOB_W,
                         FUN = object@aggregate_subsets)

    OOB_C_hat
  }

