#' Sliding similarity via Pearson's correlation
#'
#' @param Y_hat_ext Predictions from the base-learners
#' across the examples.
#'
#' @param lambda window size for computing correlations
#'
#' @return a list with a correlation matrix
#' for each prediction point
#'
#' @keywords internal
#'
#' @export
sliding_similarity <-
  function(Y_hat_ext, lambda) {
    if (is.null(colnames(Y_hat_ext))) {
      stop("null colnames in sliding similarity")
    }

    seq. <- seq_len(nrow(Y_hat_ext))[-seq_len(lambda)]

    sld_sim <- vector("list", nrow(Y_hat_ext) - lambda)
    for (i in seq.) {
      sliding_Yhat <- Y_hat_ext[(i-lambda+1):i, ]

      avg_ <- mean(unlist(sliding_Yhat), na.rm = TRUE)

      correlation_matrix <- stats::cor(as.matrix(sliding_Yhat))
      if (any(is.na(correlation_matrix))) {
        correlation_matrix[is.na(correlation_matrix)] <- 1.
      }

      sld_sim[[i-lambda]] <- correlation_matrix
    }
    sld_sim
  }


#' Sequential Re-weighting for controlling predictions' redundancy
#'
#' Besides ensemble heterogeneity we encourage diversity
#' explicitly during the aggregation of the output of experts.
#' This is achieved by taking into account not only predictions
#' of performance produced by the arbiters, but also the
#' correlation among experts in a recent window of observations.
#'
#' @param sliding_similarity list of pairwise similarity values. See
#' \code{\link{sliding_similarity}}
#'
#' @param W weights before re-weighting
#'
#' @keywords internal
#'
#' @export
sequential_reweighting <-
  function(sliding_similarity, W) {
    nws <- colnames(W)

    seq. <- seq_len(nrow(W))
    for (j in seq.) {
      W_S <- unlist(W[j, ])

      new_order <- order(W_S, decreasing = TRUE)
      first_order <- order(new_order)

      W_S <- proportion(W_S)
      names(W_S) <- nws
      W_S <- W_S[new_order]

      W_redistd <- weight_redist(sliding_similarity[[j]], W_S)
      W_redistd <- W_redistd[first_order]

      W[j, ] <- W_redistd
    }
    W
  }

weight_redist <-
  function(ssimilarity, W) {
    W_final <- rep(NA_real_, times = length(W))
    names(W_final) <- names(W)
    W_final[1] <- W[1]

    for (i in seq_along(W)[-1]) {
      W_i <- W[i]
      nW_i <- names(W_i)

      RankedW_at_i <- names(W_final[!is.na(W_final)])

      W_final[i] <- W_i
      for (exprt in RankedW_at_i) {
        cor_ij <- ssimilarity[nW_i, exprt]
        eta_ij <- cor_ij * W_final[[exprt]] * W_final[i]

        W_final[[exprt]] <- W_final[[exprt]] + eta_ij
        W_final[i] <- W_final[i] - eta_ij
      }
    }
    W_final
  }
