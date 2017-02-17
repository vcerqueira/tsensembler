setMethod("show",
          signature("tse_hat"),
          function(object) {
            null.slot <- sapply(
              sapply(slotNames(object), function(i) slot(object, i)),
              is.null
            )

            avai.slots <- paste(names(which(!null.slot)), collapse = ", ")

            cat("Predictions produced by a tseModel\n")
            cat("Available slots are:", avai.slots, "\n\n")

            names_M <- names(object@Y_hat)
            len <- length(names_M)

            unique_M <- unique(vcapply(names_M, function(j) .splitBy_(j)[1]))
            unique_M <- paste(unique_M, collapse = ", ")

            cat("The predictions were calculated with", len, "base models\n")
            cat("Base models include", unique_M, "\n")
          }
)

#' Proportionalize the predictions of base models
#'
#' @param object object of class \code{tse_hat}
#'
#' @export
setGeneric("prop_hat", function(object) {
  standardGeneric("prop_hat")
})

#' Proportionalize the predictions of base models
#'
#' @param object object of class \code{tse_hat}
#'
#' @export
setMethod("prop_hat",
          signature("tse_hat"),
          function(object) {
            Y_hat <- lapply(object@Y_hat,
                            function(y_hat) y_hat[["preds"]])
            as.data.frame(Y_hat)
          }
)


#' The individual models compute point-wise loss
#'
#' @param object object of class \code{tse_hat}, the predictions
#' made by a \code{tseModel} object.
#' @param prop Logical for proportionalization. If \strong{TRUE}  
#' results are return as \code{data.frame}. Otherwise returns
#' as \code{list}. Defaults to TRUE
#' @param lossFUN loss function to compute across base learners
#'
#' @export
setGeneric("loss_M", function(object, prop = TRUE, lossFUN = se) {
  standardGeneric("loss_M")
})

#' The individual models compute point-wise loss
#'
#' @param object object of class \code{tse_hat}, the predictions
#' made by a \code{tseModel} object.
#' @param prop Logical for proportionalization. If \strong{TRUE}  
#' results are return as \code{data.frame}. Otherwise returns
#' as \code{list}. Defaults to TRUE
#' @param lossFUN loss function to compute across base learners
#'
#' @export
setMethod("loss_M",
          signature("tse_hat"),
          function(object, prop = TRUE, lossFUN = se) {
            loss <- lapply(object@Y_hat,
                           function(x) lossFUN(x[["preds"]], 
                                               x[["trues"]]))
            if (prop) loss <- as.data.frame(loss)
            
            loss
          }
)

