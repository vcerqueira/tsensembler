#' CA generaliser using arbitrage
#'
#' @param obj a \code{\link{constructive_aggregation-class}}
#' object class
#' @param predict_info predict information -- output
#' of \code{\link{hat_info}} function
#' @param newdata newdata to predict
#'
#' @keywords internal
CA.ADE_hat <-
  function(obj, predict_info, newdata) {
    C_hat_OOB <- OOB_Chat(obj)
    Y_OOB <- get_y(obj@out_of_bag$oob, obj@form)
    OOB <- obj@out_of_bag$oob
    tgt <- get_target(obj@form)

    OOB <- subset(OOB, select = -which(colnames(OOB) %in% tgt))

    E_OOB <- base_models_loss(C_hat_OOB,Y_OOB,ae)

    train_metadata <-
      lapply(E_OOB,
             function(l_hat) {
               cbind.data.frame(OOB,
                                score = l_hat)
             })

    meta_models <-
      lapply(train_metadata,
             function(meta_set) {
               if (any(is.na(meta_set))) {
                 meta_set <- soft.completion(meta_set)
               }
               loss_meta_learn(score ~ ., meta_set, "randomforest")
             })

    E_hat <- lapply(meta_models,
                    function(o) {
                      meta_predict(o, newdata, "randomforest")
                    })

    Y <- get_y(newdata, obj@form)

    names(E_hat) <- colnames(predict_info$C_hat)
    E_hat <- abs(as.data.frame(E_hat))

    W <- t(apply(E_hat, 1, model_weighting, "linear"))

    C <- build_committee(predict_info$C_hat, Y, 50, .5)

    combine_predictions(predict_info$C_hat, W, C)
  }

#' CA generaliser using exponentially weighted average
#'
#' @param obj a \code{\link{constructive_aggregation-class}}
#' object class
#' @param predict_info predict information -- output
#' of \code{\link{hat_info}} function
#'
#' @keywords internal
#'
#' @import opera
CA.EWA_hat <-
  function(obj, predict_info) {
    C_hat <- predict_info$C_hat

    if (ncol(C_hat) < 2) {
      return(as.vector(unlist(C_hat)))
    }

    C_hat_OOB <- OOB_Chat(obj)
    Y_OOB <- get_y(obj@out_of_bag$oob, obj@form)

    Y_all <- c(Y_OOB, predict_info$Y)
    C_hat_all <- as.matrix(rbind(C_hat_OOB, C_hat))

    if (any(is.na(C_hat_all))) {
      C_hat_all <- soft.completion(C_hat_all)
    }

    EWA_mix <- mixture(model = "EWA", loss.type = "absolute")
    for (i in 1:length(Y_all)) {
      EWA_mix <- predict(EWA_mix, newexperts = C_hat_all[i, ], newY = Y_all[i])
    }

    W <- EWA_mix$weights[-seq_len(nrow(C_hat_OOB)),]

    if (is.null(ncol(W))) {
      W <- t(t(W))
    }


    combine_predictions(C_hat, W, NULL)
  }

#' CA generaliser using OGD
#'
#' @inheritParams CA.EWA_hat
#'
#' @keywords internal
#'
#' @import opera
CA.OGD_hat <-
  function(obj, predict_info) {
    C_hat <- predict_info$C_hat

    if (ncol(C_hat) < 2) {
      return(as.vector(unlist(C_hat)))
    }


    C_hat_OOB <- OOB_Chat(obj)
    Y_OOB <- get_y(obj@out_of_bag$oob, obj@form)

    Y_all <- c(Y_OOB, predict_info$Y)
    C_hat_all <- as.matrix(rbind(C_hat_OOB, C_hat))

    if (any(is.na(C_hat_all))) {
      C_hat_all <- soft.completion(C_hat_all)
    }

    OGD_mix <- mixture(model = "OGD", loss.type = "absolute")
    for (i in 1:length(Y_all)) {
      OGD_mix <- predict(OGD_mix, newexperts = C_hat_all[i, ], newY = Y_all[i])
    }

    W <- OGD_mix$weights[-seq_len(nrow(C_hat_OOB)),]
    if (is.null(ncol(W))) {
      W <- t(t(W))
    }

    combine_predictions(C_hat, W, NULL)
  }

#' CA generaliser using fixed share
#'
#' @inheritParams CA.EWA_hat
#'
#' @keywords internal
#'
#' @import opera
CA.FixedShare_hat <-
  function(obj, predict_info) {
    C_hat <- predict_info$C_hat

    if (ncol(C_hat) < 2) {
      return(as.vector(unlist(C_hat)))
    }


    C_hat_OOB <- OOB_Chat(obj)
    Y_OOB <- get_y(obj@out_of_bag$oob, obj@form)

    Y_all <- c(Y_OOB, predict_info$Y)
    C_hat_all <- as.matrix(rbind(C_hat_OOB, C_hat))

    if (any(is.na(C_hat_all))) {
      C_hat_all <- soft.completion(C_hat_all)
    }

    FS_mix <- mixture(model = "FS", loss.type = "absolute")
    for (i in 1:length(Y_all)) {
      FS_mix <- predict(FS_mix,
                        newexperts = C_hat_all[i, ],
                        newY = Y_all[i])
    }

    W <- FS_mix$weights[-seq_len(nrow(C_hat_OOB)),]

    if (is.null(ncol(W))) {
      W <- t(t(W))
    }

    combine_predictions(C_hat, W, NULL)
  }

#' CA generaliser using ridge regression
#'
#' @inheritParams CA.EWA_hat
#'
#' @keywords internal
#'
#' @import opera
CA.Ridge_hat <-
  function(obj, predict_info) {

    C_hat <- predict_info$C_hat

    if (ncol(C_hat) < 2) {
      return(as.vector(unlist(C_hat)))
    }

    C_hat_OOB <- OOB_Chat(obj)
    Y_OOB <- get_y(obj@out_of_bag$oob, obj@form)

    Y_all <- c(Y_OOB, predict_info$Y)
    C_hat_all <- as.matrix(rbind(C_hat_OOB, C_hat))

    if (any(is.na(C_hat_all))) {
      C_hat_all <- soft.completion(C_hat_all)
    }

    Ridge_mix <- mixture(model = "Ridge", loss.type = "square")
    Ridge_mix <-
      predict(Ridge_mix,
              newexpert = C_hat_all,
              newY = Y_all,
              online = TRUE)

    W <- Ridge_mix$weights[-seq_len(nrow(C_hat_OOB)),]

    if (is.null(ncol(W))) {
      W <- t(t(W))
    }

    combine_predictions(C_hat, W, NULL)
  }

#' CA generaliser using polynomial weighted average
#'
#' @inheritParams CA.EWA_hat
#'
#' @keywords internal
#'
#' @import opera
CA.MLpol_hat <-
  function(obj, predict_info) {
    C_hat <- predict_info$C_hat

    if (ncol(C_hat) < 2) {
      return(as.vector(unlist(C_hat)))
    }


    C_hat_OOB <- OOB_Chat(obj)
    Y_OOB <- get_y(obj@out_of_bag$oob, obj@form)

    Y_all <- c(Y_OOB, predict_info$Y)
    C_hat_all <- as.matrix(rbind(C_hat_OOB, C_hat))

    if (any(is.na(C_hat_all))) {
      C_hat_all <- soft.completion(C_hat_all)
    }

    MLpol0 <- mixture(model = "MLpol", loss.type = "absolute")
    for (i in 1:length(Y_all)) {
      MLpol0 <- predict(MLpol0, newexperts = C_hat_all[i, ], newY = Y_all[i])
    }

    W <- MLpol0$weights[-seq_len(nrow(C_hat_OOB)),]

    if (is.null(ncol(W))) {
      W <- t(t(W))
    }

    combine_predictions(C_hat, W, NULL)
  }
