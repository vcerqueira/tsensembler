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
embed.timeseries <- function(timeseries, embedding.dimension) {
  if (missing(timeseries)) stop("Provide a time series.", call. = FALSE)

  if (!is.xts(timeseries) && class(timeseries) == "numeric") {
    timeseries <- xts(timeseries,
                      order.by = seq.Date(from = Sys.Date(),
                                          by = 1,
                                          length.out = length(timeseries)))
  }
  ts.index  <- index(timeseries)[-(1:(embedding.dimension-1))]

  ts.embed  <- embed(timeseries, dimension = embedding.dimension)
  colnames(ts.embed) <- c('target', paste0('Tm', 1:(embedding.dimension-1)))

  df <- as.data.frame(xts(ts.embed, order.by = ts.index))

  df
}

#' Computing Summary Statistics on the Embedding Vectors
#'
#' This function computes summary statistics on the embedding vectors
#' The objetive is to augment the information on the recent dynamics of
#' the series.
#'
#' @param data a matrix or data.frame composed by the embedding vectors
#' @param stats.only logical. if TRUE, it returns the stats as data.frame. Otherwise
#' stats are returned cbinded with the embedded series. Defaults to FALSE.
#'
#' @seealso \code{\link{embed.timeseries}} for details on the embedding
#' of a time series into an Euclidean space.
#'
#' @examples
#' \dontrun{
#' require(xts)
#' ts <- as.xts(rnorm(100L), order.by = Sys.Date() + rnorm(100L))
#' embedded.ts <- embedStats(embed.timeseries(ts, 20L))
#' }
#'
#' @return A data.frame with embedding vectors binded with summary statistics
#'
#' @export
embedStats <- function(data, stats.only = FALSE) {
  if (methods::is(data, "vector")) return(c(data, mean = mean(data), sd = sd(data)))

  dstats <- vapply(c(mean, sd), function(stat) apply(data, 1, stat), double(NROW(data)))
  colnames(dstats) <- c("mean", "sd")

  if (stats.only)
    res <- dstats
  else
    res <- cbind.data.frame(data, dstats)

  res
}

#' unembed.timeseries
#'
#' This function reverts the process of embed.timeseries function.
#' It takes an embedded series and returns the series as a sequential vector.
#'
#' @param train train series data.frame
#' @param test test series data.frame
#'
#' @return a time series
#' @export
unembed.timeseries <- function(train, test = NULL) {
  train <- unlistn(c(as.vector(train[1, NCOL(train):1]), as.vector(train[-1,1])))

  if (!is.null(test)) {
    test  <- unlistn(as.vector(test[, 1]))
    data  <- c(train, test)
  } else
    data <- train

  list(data = data, train = train, test = test)
}

HURST <- function(x) {
  cwtwnoise <- Rwave::DOG(x, 10, 3, 1, plot=FALSE)
  mcwtwnoise <- Mod(cwtwnoise)
  mcwtwnoise <- mcwtwnoise * mcwtwnoise
  wspwnoise <- Rwave::tfmean(mcwtwnoise, plot=FALSE)

  Rwave::hurst.est(wspwnoise, 1:7, 3, plot = FALSE)[[2]]
}

#' Time series dynamics
#'
#' @param data Numeric matrix to compute dynamics in, such as an
#' embedded time series. 
#'
#' @export
ts.dynamics <- function(data) {
  seq. <- seq_len(nrow(data))
  seq_no1 <- seq.[-1]
  K <- ncol(data)

  trd <- apply(data, 1, trend)
  skw <- apply(data, 1, moments::skewness)
  kts <- apply(data, 1, moments::kurtosis)
  mle <- apply(data, 1, function(r) {
    Reduce(max,
           nonlinearTseries::divergence(
             nonlinearTseries::maxLyapunov(time.series = r,
                                           min.embedding.dim = ceiling(K / 4),
                                           max.embedding.dim = ceiling(K / 2),
                                           radius = ceiling(K / 6),
                                           do.plot = FALSE)
           )
    )
  })

  hrst <- apply(data, 1, HURST)
  selfs <- apply(data, 1, function(j) tseries::terasvirta.test(x = as.ts(unlistn(j)), lag = 3)$p.value)
  serialcorr <- apply(data, 1, function(j) Box.test(j)$p.val)

  dStats <- data.frame(trd, skw, kts, mle, hrst, selfs, serialcorr)
  dStats[is.na(dStats)] <- double(1L)
  colnames(dStats) <- c("trd", "skw", "kts", "mle", "hrst", "selfs", "serialcorr")

  nzv_cols <- caret::nearZeroVar(dStats)
  if (length(nzv_cols) > 0L) {
    dStats <- subset(dStats, select = -nzv_cols)
  }
  rownames(dStats) <- NULL

  preproc <- caret::preProcess(dStats)
  dStats <- predict(preproc, dStats)

  dStats
}

#' Normalized variance by range
#'
#' @param x time series as vector
#' @param ... further params to \code{var}, \code{max} and \code{min}
#'
#' @export
rangedvar <- function(x, ...) var(x, ...) / (max(x, ...) - min(x, ...))

#' trend
#'
#' @param x x
#' @param ... further params to \code{sd}
#'
#' @export
trend <- function(x, ...) sd(x, ...) / sd(diff(x), ...)
