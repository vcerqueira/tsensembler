#' pointwise metalearn with rpart
#'
#' @param metadata metadata from \code{setup.metadata}
#' @param timeseq time sequence of test set
#' @param steps metamodels are retrained every \code{steps}
#'
#' @import rpart
#' 
#' @export
point.metalearn.rpart <- function(metadata, timeseq, steps = 10) {
  len <- length(timeseq)
  oM_hat <- lapply(metadata, function(metaset) {
    retrain.time <- seq(from = timeseq[1], to = timeseq[len], by = steps)
    om_yhat <- numeric(len)
    for (j in seq_along(timeseq)) {
      N_j <- timeseq[j]
      if (N_j %in% retrain.time) {
        capture.output(metamodel <- rpart(score ~., metaset[seq_len(N_j), ]))
      }
      om_yhat[j] <- predict(metamodel, metaset[N_j + 1L, ])
    }
    om_yhat
  })
  as.data.frame(oM_hat)
}

#' point metalearn with GP
#'
#' @inheritParams point.metalearn.rpart
#'
#' @import kernlab
#'
#' @export
point.metalearn.gausspr <- function(metadata, timeseq, steps = 10) {
  len <- length(timeseq)
  oM_hat <- lapply(metadata, function(metaset) {
    retrain.time <- seq(from = timeseq[1], to = timeseq[len], by = steps)
    om_yhat <- numeric(len)
    for (j in seq_along(timeseq)) {
      N_j <- timeseq[j]
      if (N_j %in% retrain.time) {
        capture.output(metamodel <- gausspr(score ~., metaset[seq_len(N_j), ], kernel = "vanilladot"))
      }
      om_yhat[j] <- predict(metamodel, metaset[N_j + 1L, ])
    }
    om_yhat
  })
  as.data.frame(oM_hat)
}

#' pointwise metalearn with random forest
#'
#' @param metadata metadata from \code{setup.metadata}
#' @param timeseq time sequence of test set
#' @param steps metamodels are retrained every \code{steps}
#'
#' @import ranger
#' @export
point.metalearn.rf <- function(metadata, timeseq, steps = 10) {
  len <- length(timeseq)
  oM_hat <- lapply(metadata, function(metaset) {
    retrain.time <- seq(from = timeseq[1], to = timeseq[len], by = steps)
    om_yhat <- numeric(len)
    for (j in seq_along(timeseq)) {
      N_j <- timeseq[j]
      if (N_j %in% retrain.time) {
        capture.output(metamodel <- ranger(score ~., metaset[seq_len(N_j), ], num.trees = 1000, write.forest=TRUE, importance = 'impurity'))
      }
      om_yhat[j] <- predict(metamodel, metaset[N_j + 1L, ])$predictions
    }
    om_yhat
  })
  as.data.frame(oM_hat)
}

#' pointwise metalearn with random forest
#'
#' @param metadata metadata from \code{setup.metadata}
#' @param timeseq time sequence of test set
#' @param steps metamodels are retrained every \code{steps}
#'
#' @import ranger
#' @export
point.metalearn.rf.augmented <- function(metadata, timeseq, steps = 10) {
  len <- length(timeseq)
  oM_hat <- lapply(metadata, function(metaset) {
    retrain.time <- seq(from = timeseq[1], to = timeseq[len], by = steps)
    om_yhat <- numeric(len)
    for (j in seq_along(timeseq)) {
      N_j <- timeseq[j]
      if (N_j %in% retrain.time) {
        capture.output(metamodel <- ranger(score ~., metaset[seq_len(N_j), ], num.trees = 1000, write.forest=TRUE, importance = 'impurity'))
      }
      om_yhat[j] <- predict(metamodel, metaset[N_j + 1L, ])$predictions
    }
    list(om_yhat = om_yhat, metamodel = metamodel)
  })
  E_hat <- lapply(oM_hat, function(j) j[[1]])
  oM <- lapply(oM_hat, function(j) j[[2]])

  list(E_hat = as.data.frame(E_hat), oM = oM)
}

#' Pointwise train of models
#'
#' @param train train set
#' @param test test set
#' @param form formula
#' @param learner learning models
#' @param learner.pars learning model parameters
#'
#' @export
pointwise_M_hat_seq <- function(train, test, form, learner, learner.pars) {
  K <- get_embedsize(train)
  target <- get_target(form)

  tr_len <- nrow(train)
  ts_len <- nrow(test)

  train_data <- rbind.data.frame(train, test)

  seq_ids <- seq_len(ts_len) - 1L
  Y_hat <- lapply(seq_ids, function(hat_point) {
    seq_train <- seq_len(tr_len + hat_point)

    seen_data <- train_data[seq_train, ]
    unseen_data <- train_data[-seq_train, ]

    cat(tr_len, "+", hat_point, "\n")
    cat(nrow(seen_data), "::", nrow(unseen_data),":::", nrow(train_data),"\n")

    M <- learnM(form, seen_data, learner, learner.pars, K)

    Y_hat <- predict(M, unseen_data)
    Y_hat.prop <- prop_hat(Y_hat)[1, ]
    Y_hat.prop
  })
  cat("Finito\n")
  
  Y_hat <- rbind_(Y_hat)
  Y_hat
}