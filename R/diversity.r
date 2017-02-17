#' Diversity analysis between pairs of models
#'
#' Methods: correlation, covariance, mutual info,
#' disagreement measure and chi-square
#'
#' ref: Measuring Diversity in Regression Ensembles
#' author of paper: haimonti duta
#'
#' @param Y_hat matrix with colwise predictions
#' @param Y vector of true values
#'
#' @export
diversity <- function(Y_hat, Y = NULL) {
  # Correlation (inversely prop)
  corY <- cor(Y_hat)
  # Covariance (inversely prop)
  covY <- cov(Y_hat)
  # Mutual Information (inversely prop)
  mi <- -.5 * log(1 - corY ^ 2)
  mi[is.infinite(mi)] <- 1.

  # Chi-square (directly prop)
  combine_cols <- as.data.frame(combn(ncol(Y_hat), 2))
  cnames <- colnames(Y_hat)

  chisq <- rbind_(
    lapply(combine_cols, function(ypair) {

      suppressWarnings(tmp <- tryCatch(
        chisq.test(Y_hat[ ,ypair[1]], Y_hat[ ,ypair[2]]), error = function(e) list()))

      c(m1 = cnames[ypair[1]],
        m2 = cnames[ypair[2]],
        pval = tmp$p.value,
        df = tmp$parameter,
        chiq = tmp$statistic)
    })
  )
  rownames(chisq) <- NULL

  adjMat_chisq <- data.frame(
    igraph::get.adjacency(
      igraph::graph.data.frame(
        chisq, directed = FALSE), attr = "pval", sparse = FALSE),
    stringsAsFactors = FALSE
  )
  adjMat_chisq[] <- lapply(adjMat_chisq, type.convert, na.strings = "")
  adjMat_chisq[is.na(adjMat_chisq)] <- 0.

  # Disagreement
  if (!is.null(Y)) {
    beta <- apply(Y_hat, 1, sd, na.rm = TRUE)
    alpha <- Y

    seq. <- seq_along(alpha)

    Y_hat.acc <- as.data.frame(
      lapply(as.data.frame(Y_hat), function(y_hat) {
        vlapply(seq., function(i) {
          y_hat[i] < alpha[i] + beta[i] && y_hat[i] > alpha[i] - beta[i]
        })
      })
    )

    dis <- rbind_(
      lapply(combine_cols, function(ypair) {
        y_1 <- Y_hat.acc[ ,ypair[1]]
        y_2 <- Y_hat.acc[ ,ypair[2]]

        dis <- (sum(!y_1 & y_2) + sum(y_1 & !y_2)) / length(y_1)

        c(m1 = cnames[ypair[1]],
          m2 = cnames[ypair[2]],
          dis = dis)
      })
    )
    rownames(dis) <- NULL

    adjMat_dis <- data.frame(
      igraph::get.adjacency(
        igraph::graph.data.frame(
          dis, directed = FALSE), attr = "dis", sparse = FALSE),
      stringsAsFactors = FALSE
    )
    adjMat_dis[] <- lapply(adjMat_dis, type.convert, na.strings = "")
    adjMat_dis[is.na(adjMat_dis)] <- 0.
  } else
    adjMat_dis <- NULL

  rm.null(
    list(correlation = corY,
         covariance = covY,
         mutual_information = mi,
         disagreement = adjMat_dis,
         chi_square = adjMat_chisq)
  )
}

#' Relative accuracy measure for regression
#'
#' This function computes the accuracy of regression models
#' according to the predictions of all models and the actual value.
#'
#' @param Y_hat predictions as df
#' @param Y true values
#'
#' @return Logical \code{data.frame} with the accuracy of each model
#'
#' @export
model_accuracy <- function(Y_hat, Y) {
  beta <- apply(Y_hat, 1, sd, na.rm = TRUE)
  alpha <- Y
  seq. <- seq_along(alpha)

  as.data.frame(
    lapply(Y_hat, function(y_hat) {
      vlapply(seq., function(i) {
        y_hat[i] < alpha[i] + beta[i] && y_hat[i] > alpha[i] - beta[i]
      })
    })
  )
}
