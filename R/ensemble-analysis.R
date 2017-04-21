#' pairedComparisons on Multiple Time Series \code{ComparisonResults} objects
#'
#' This function applies the \code{\link[performanceEstimation]{pairedComparisons}} method
#' to a list of \code{ComparisonResults} objects. This can be used to analyze differences
#' in performance of two or more models accross several predictive tasks.
#'
#' @param estimationResults A list of \code{ComparisonResults} objects.
#' @param baseline Character name of workflow to compared methods to.
#'
#' @seealso Check \code{\link[performanceEstimation]{pairedComparisons}} for a comprehensive
#' explanation of the core method.
#'
#' @return Paired Comparison results accross the several predictive tasks.
#' @export
Multi.pairedComparisons <- function(estimationResults, baseline) {
  pC <- lapply(estimationResults, function(r) {
    performanceEstimation::pairedComparisons(obj = r, baseline = baseline)
  })
  names(pC) <- paste0('Experiment_', seq_along(estimationResults))
  pC
}

#' Statistical Wins and Losses
#'
#' WinLoss is an utility function built of top of \code{pairedComparisons}
#' method (\code{performanceEstimation} package) that analyses the
#' performance of models. Particularly, \code{WinLoss} probes for statistically
#' significant differences in the performance of models.
#'
#' @param pairedComp Paired Comparison results. This is a list
#' with experiments computed from the mean and median scores of
#' workflows.
#' @param metric Metric used to compare methods. Defaults to Mean Squared Error (MSE).
#' @param test Statistical significance test to perform. Defaults to \code{WilcoxonSignedRank.test}.
#' The other option is the \code{t.test}
#' @param pretty Logical: If TRUE, results are rendered in a more visually appealling way.
#' @param pvalue Significance threshold
#'
#' @seealso \code{\link[performanceEstimation]{pairedComparisons}}
#'
#' @examples
#'
#' \dontrun{
#' pairedComps <- pairedComparisons(ComparisonResults.obj, baseline)
#' wl <- WinLoss(pairedComps, "WilcoxonSignedRank.test")
#' }
#'
#' @return a data.frame with the significance test results accross several experiments
#'
#' @export
WinLoss <- function(pairedComp,
					          metric = "mse",
					          test = "t.test",
                    pretty = FALSE,
                    pvalue = 0.05) {
  if (!test %in% c("WilcoxonSignedRank.test", "t.test")) {
    stop("Please select a valid test (either \"t.test\" or \"WilcoxonSignedRank.test\")")
  }
  pC <- pairedComp[[metric]]
  baseline <- pC[["baseline"]]
  WSR <- as.data.frame(pC[[test]])
  methods <- row.names(WSR)
  DT <- data.frame(methods, WSR)
  if (test == "WilcoxonSignedRank.test") colnames(DT) <- c("Method", "MedScore", "DiffMedScore", "pval")
  if (test == "t.test") colnames(DT) <- c("Method", "AvgScore", "DiffAvgScores", "pval")
  rownames(DT) <- NULL

  methods <- DT[,"Method"]
  len <- length(methods)

  WL <- data.frame(Win    = rep(NA, len),
                   SigWin = rep(NA, len),
                   Loss   = rep(NA, len),
                   SigLoss= rep(NA, len))
  rownames(WL) <- methods
  for (i in seq_along(methods)) {
    if (methods[i] == baseline) next
    WL[i, ] <- as.integer(c(DT[i, 3] > 0,
                            DT[i, 3] > 0 && DT[i, "pval"] < pvalue,
                            DT[i, 3] < 0,
                            DT[i, 3] < 0 && DT[i, "pval"] < pvalue))
  }
  if (pretty) WL <-prettify(WL)
  WL
}

#' Statistical Wins and Losses on Multiple Time Series
#'
#' Utility function built on top of \code{WinLoss} function. This function applies and combines
#' the results of \code{WinLoss} accross several experiments.
#'
#' @param ResultsList Experiment results, process by \code{Multi.pairedComparisons} function
#' @param test Statistical significance test to perform. Defaults to \code{WilcoxonSignedRank.test}.
#' The other option is the \code{t.test}
#' @param pretty Logical: If TRUE, results are rendered in a more visually appealling way.
#' @param ... Further arguments to pass to \code{WinLoss} function.
#'
#' @seealso \code{\link{Multi.pairedComparisons}}
#'
#' @examples
#'
#' \dontrun{
#'
#' pairedComps <- Multi.pairedComparisons(list_of_ComparisonResults, baseline)
#' wl <- Multi.WinLoss(pairedComps, "WilcoxonSignedRank.test")
#'
#' }
#'
#' @return a data.frame with significance results.
#'
#' @export
Multi.WinLoss <- function(ResultsList, test = "t.test", pretty = TRUE, ...) {
  WL <- Reduce('+',
               lapply(ResultsList, function(result) {
                 WinLoss(result, test = test, pretty = FALSE, ...)
                 })
               )
  if (pretty) WL <- prettify(WL)

  WL
}

#' prettify
#'
#' Utility function that transforms the results from \code{WinLoss}
#' into a more visual appealing frame, to facilitate analysis.
#' @param WL Results from Significance experiments
#'
#' @return Pretty results from Significance experiments
prettify <- function(WL) {
  strWL <- data.frame(rep(NA, NROW(WL)))
  colnames(strWL) <- rownames(WL)[which(is.na(WL$Win))]

  for (i in seq_len(NROW(WL))) {
    if (length(which(is.na(WL[i, ]))) == 0 )
      strWL[i, ] <- paste0(WL[i, "Win"]," (",
                           WL[i, "SigWin"],") / ",
                           WL[i, "Loss"]," (",
                           WL[i, "SigLoss"],")")
  }
  rownames(strWL) <- rownames(WL)

  na.omit(strWL)
}

#' Average Rank of Workflows
#'
#' This function takes as input a paired comparisons object
#' of workflows across several experiments and computes the
#' average rank and respective deviation of each workflow.
#'
#' @param obj Paired Comparisons results
#' @param metric Metric used to rank workflows. Defaults to \strong{rmse}
#'
#' @return Data.frame with Average and Standard deviation Rank statistics
#' of workflows
#'
#' @export
avg_ranks <- function(obj, metric = "rmse") {
  ranks <- as.data.frame(lapply(obj,
                                function(z) {
                                  z[[metric]][["avgRksWFs"]]
                                  }))

  data.frame(workflow = rownames(ranks),
             meanRank = round(rowMeans(ranks), 2),
             sdRank = round(apply(ranks, 1, sd), 2), row.names = NULL)
}


#' worflows err over time
#'
#' @param compRes ComparisonResults class object
#'
#' @export
wf_err <- function(compRes) {
  sapply(compRes, function(wf) {
    rowMeans(
      sapply(wf@iterationsInfo, function(x) {
        ae(x[["trues"]], x[["preds"]])
      }), na.rm = TRUE)
  })
}

