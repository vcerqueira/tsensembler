setGeneric("getTSE", function(object, .iter = 1L) {
  standardGeneric("getTSE")
})

#' Getting the best TSE model from a given Monte Carlo iteration
#'
#' This utility function handles a \code{ComparisonResults} object to
#' extract the best performing learning ensemble model in a given \code{.iter}
#' Monte Carlo repetition.
#'
#' @param object \code{ComparisonResults} object
#' @param .iter Monte Carlo repetition
#'
#' @return A tseModel-class object
setMethod("getTSE",
          signature("ComparisonResults"),
          function(object, .iter) {
            .bestWF <- topPerformer(object, "rmse", names(object))
            .bestID <- as.numeric(which(sapply(object[[1]], function(x) {
              x@workflow@name
            }) == .bestWF@name))

            .gII <- getIterationsInfo(object, workflow = .bestID)
            .l <- length(.gII)
            if (.iter > .l) {
              stop("Invalid iteration. How many Monte Carlo reps are run?")
            }
            .gII[[.iter]]$tse
          }
)