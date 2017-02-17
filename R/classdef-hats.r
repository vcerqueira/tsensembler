setClassUnion("OptionalNumeric", c("numeric","NULL"))
setClassUnion("OptionalList", c("list","NULL"))
setClassUnion("OptionalDF", c("data.frame","NULL"))

#' Predictions produced by a \code{tseModel} object
#'
#' @slot y_hat Combined predictions. Numeric vector
#' @slot Y_hat List of predictions made by individual models
#' @slot Y_W Normalized weights of individual models
#' @slot Y_committee List of top models at each prediction step
#' @slot Y_rank List of rank of models at each prediction step
#' @slot data data used to produce predictions
#'
#' @seealso \code{\link{tseModel-class}} for the \code{tseModel} class.
#'
#' @return Set of results from predicting with a \code{tseModel}
#' ensemble
#'
#' @export
setClass("tse_hat",
         slots = c(y_hat = "numeric",
                   Y_hat = "OptionalList",
                   Y_W = "OptionalDF",
                   Y_committee = "OptionalList",
                   Y_rank = "OptionalList",
                   data = "OptionalDF")
)

#' Predictions by a tseModel
#'
#' @param y_hat Combined predictions. Numeric vector
#' @param Y_hat List of predictions made by individual models
#' @param Y_W Normalized weights of individual models
#' @param Y_committee List of top models at each prediction step
#' @param Y_rank List of rank of models at each prediction step
#' @param data data used to produce predictions
#'
#' @seealso \code{\link{tseModel-class}} for the \code{tseModel} class.
#'
#' @return Set of results from predicting with a \code{tseModel}
#' ensemble
#'
#' @export
tse_hat <- function(y_hat = "numeric",
                    Y_hat = "OptionalList",
                    Y_W = "OptionalDF",
                    Y_committee = "OptionalList",
                    Y_rank = "OptionalList",
                    data = "data.frame") {

  if (missing(Y_hat)) Y_hat <- NULL
  if (missing(Y_rank)) Y_rank <- NULL
  if (missing(Y_committee)) Y_committee <- NULL
  if (missing(data)) data <- NULL
  if (missing(Y_W)) Y_W <- NULL

  new("tse_hat",
      y_hat = y_hat,
      Y_hat = Y_hat,
      Y_W = Y_W,
      Y_committee = Y_committee,
      Y_rank = Y_rank,
      data = data)
}

