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
        capture.output(metamodel <- gausspr(score ~., metaset[seq_len(N_j), ]))
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
        capture.output(metamodel <- ranger(score ~., metaset[seq_len(N_j), ], num.trees = 750, write.forest=TRUE))
      }
      om_yhat[j] <- predict(metamodel, metaset[N_j + 1L, ])$predictions
    }
    om_yhat
  })
  as.data.frame(oM_hat)
}
