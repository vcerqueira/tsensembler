#' Contiguity check
#'
#' This function checks the contiguity of
#' out-performance of predictive models,
#' and respective regions where this occurs
#'
#' @param C committee set
#' @param alpha contiguity interval size
#'
#' @keywords internal
#'
#' @export
contiguous_count <-
  function(C, alpha) {
    C_sign <-
      Map(function(x) paste(x, collapse = "_"), C)

    C_sign <- unlist(C_sign)

    C_rle <- rleid(C_sign)
    ids <- which(table(C_rle) >= alpha)

    covered_regions <-
      lapply(ids,
             function(x) {
               which(C_rle %in% x)
             })

    names_regions <-
      sapply(ids, function(x) {
        C_sign[min(which(x == C_rle))]
      })

    names(covered_regions) <- names_regions

    top_subsets <- as.list(
      unique(
        C_sign[C_rle %in% ids]
      )
    )

    list(top_subsets = top_subsets,
         covered_regions = covered_regions)
  }

#' Prune subsets by out-performance
#'
#' @param loss_info loss data in unseen observations
#' @param lambda smoothing window size
#' @param depth depth size how large is the maximum size
#' of the subsets. If NULL, defaults to no. of predictors
#' minus one.
#'
#' @keywords internal
#'
#'
#' @export
prune_c_outperformance <-
  function(loss_info, lambda, depth = NULL) {
    rolled_loss <- roll_mean_matrix(loss_info, lambda)

    N <- ncol(loss_info)

    if (is.null(depth)) {
      seq_N <- 1:(N-1)
    } else {
      seq_N <- depth
    }

    C_all <-
      lapply(seq_N,
             function(i) {
               l1apply(rolled_loss,
                       function(j) {
                         point_loss <- unlist(rolled_loss[j, ])
                         order(point_loss)[seq_len(i)]
                       })
             })

    C_all <- lapply(C_all, rm.len0)

    rm.len0(C_all)
  }


#' Prune subsets by contiguity
#'
#' @param C list of outperformers. The output of
#' \code{\link{prune_c_outperformance}}
#' @param alpha threshold for cut-off. integer
#'
#' @seealso \code{\link{prune_c_outperformance}},
#' whose output is the input for this function
#'
#' @keywords internal
#'
#' @export
prune_c_contiguity <-
  function(C, alpha) {
    C <- lapply(C, contiguous_count, alpha)
    covered_regions <- lapply(C, function(x) x[[2]])
    C <- lapply(C, function(x) x[[1]])

    C <- unlist(C, recursive = FALSE)

    C <-
      lapply(C,
             function(x) {
               as.integer(split_by_(x))
             })

    if (length(C) < 1) {
      warning("cant find C's with these pars.")
    } else if (length(C) == 1) {
      warning("a single C found.")
    } else {}

    list(C = rm.len0(C),
         covered_regions = covered_regions)
  }

