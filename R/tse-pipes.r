#' Get true values from test set
#' from formula and the test set
#' 
#' @param test test set
#' @param form formula
#' 
#' @export
get_y <- function(test, form) model.response(model.frame(form, test, na.action = NULL))

#' get cols that are part of the embedding vector
#'
#' @param x embedded time series
#'
#' @export
get_embedcols <- function(x) grepl("^Tm[0-9]?[0-9]$", colnames(x))

#' Get embedding dimension from an embedded time series
#'
#' @param x embedded time series
#'
#' @export
get_embedsize <- function(x) length(which(get_embedcols(x))) + 1L

#' readSignal
#'
#' readSignal is an utility function used to munge a new data point into
#' the \code{data.frame} format identic to a given embedded time series.
#' Using this application turns a \code{raw_signal} into a format suitable for
#' tseModel objects to operate.
#'
#' @param raw_signal a new embedding vector to process into prediction-ready format
#'
#' @return A munged embedding vector
readSignal <- function(raw_signal) {
  names(raw_signal) <- paste0("V", 2:(length(raw_signal)+1L))
  #sig. <- t(embedStats(raw_signal))

  cbind.data.frame(target = -1L, as.data.frame(t(raw_signal)))
}

#' Assembling Workflows Together
#'
#' This is an utility function used to assemble workflows by embedding size.
#' The goal is the to apply paired comparisons procedures to those workflows
#'
#' @param ensembleResults experiment results
#'
#' @return grouped experiment results
.assemble <- function(ensembleResults) {
  res <- list()
  res[[1]] <- ensembleResults[[1]]
  seq. <- seq_along(ensembleResults)

  for (z in seq.[-1]) {
    res[[1]][[1]] <- c(res[[1]][[1]], ensembleResults[[z]][[1]])
  }

  names(res[[1]][[1]]) <- paste(names(res[[1]][[1]]), seq., sep = "_")

  res[[1]]
}
