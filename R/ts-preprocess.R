#' Embedding a Time Series
#'
#' This function embeds a time series into an Euclidean space.
#' This implementation is based on the function \code{embed} of
#' \strong{stats} package and has theoretical backgroung on
#' reconstruction of attractors (see Takens, 1981).
#' This shape transformation of the series allows for
#' the use of any regression tool available to learn
#' the time series. The assumption is that there are no long-term
#' dependencies in the data.
#'
#' @param timeseries a time series of class \"xts\".
#' @param embedding.dimension an integer specifying the embedding dimension.
#'
#' @return An embedded time series
#'
#' @seealso \code{\link[stats]{embed}} for the details of the embedding procedure.
#'
#' @examples
#' \dontrun{
#' require(xts)
#' ts <- as.xts(rnorm(100L), order.by = Sys.Date() + rnorm(100L))
#' embedded.ts <- embed.timeseries(ts, 20L)
#' }
#'
#' @export
#'
#' @import zoo
#' @import xts
embed_timeseries <- function(timeseries, embedding.dimension) {
  if (missing(timeseries)) stop("Provide a time series.", call. = FALSE)

  if (!is.xts(timeseries) && class(timeseries) == "numeric") {
    timeseries <- xts(timeseries,
                      order.by = seq.Date(from = Sys.Date(),
                                          by = 1,
                                          length.out = length(timeseries)))
  }
  ts.index  <- index(timeseries)[-(1:(embedding.dimension-1))]

  ts.embed  <- stats::embed(timeseries, dimension = embedding.dimension)
  colnames(ts.embed) <- c('target', paste0('Tm', 1:(embedding.dimension-1)))

  df <- as.data.frame(xts(ts.embed, order.by = ts.index))

  df
}

#' Soft Imputation
#' 
#' @param x data
#' 
#' @import softImpute
#' 
#' @export
soft.completion <- function(x) {
  if ("data.frame" %in% class(x)) x <- as.matrix(x)
  as.data.frame(complete(x, softImpute(x)))
}