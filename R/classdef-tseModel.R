setOldClass("xts")
setOldClass("zoo")
setClassUnion("tsSource", c("list", "xts", "zoo"))
setClassUnion("OptionalNumeric", c("numeric","NULL"))
setClassUnion("OptionalCharacter", c("character","NULL"))
setClassUnion("OptionalLogical", c("logical","NULL"))
setClassUnion("OptionalList", c("list","NULL"))

#' tseModel-class
#'
#' tseModel is an S4 class that contains the ensemble model. Besides
#' the base learning algorithms -- \code{baseModels} -- tseModel class
#' contains information about other meta-data used to compute predictions
#' for new upcoming data.
#'
#' @slot baseModels list comprising the base learners.
#'
#' @slot preWeights Normalized relative weights of the base learners according to
#' their performance on the available data.
#'
#' @slot N Number of different base Models.
#'
#' @slot modelDist base learner distribution with respect to the type of learner.
#' That is, the number of Decision Trees, SVMs, etc.
#'
#' @slot form Formula of class formula. Essentially used for reference.
#'
#' @slot ma.N The number of periods to average over the Squared Error when
#' computing \strong{MASE} (\emph{Moving Average Squared Error}). This parameter is only
#' used when the function applied to weight the base models is based on
#' MASE. The rationale behind this simple heuristic is that base models are weighted according
#' to their recent performance. Here, recent performance is formalized and quantified as
#' inversely proportional to the models' MASE. The dynamics of the moving average
#' yield a flexibility to the combined model, in the sense that it is self-adaptable when
#' concept drift occurs. Particularly, this parameter \code{ma.N} controls
#' the reactiveness of the learning system to such events.
#' A small value of P leads to greater reactiveness, but also makes the ensemble
#' susceptible to be fooled by outliers. Conversely, higher values of \code{ma.N}
#' lead to greater stability, while sacrificing some responsiveness. This trade-off
#' is known in the literature as the stability-plasticity dilemma (Carpenter et al., 1991).
#'
#' @slot embedding.dimension The maximum embedding dimension used to transform the series.
#'
#' @slot committee.ratio A numeric value between 0 and 1 representing the ratio
#' of base learners that comprise the committee at each prediction time.
#' If \code{committee.ratio} equals, say, 0.2, at time \emph{t}, the ensemble will the
#' 20\% best base learners up to time \emph{t - 1}.
#'
#' @slot aggregationFUN The function name used to combine the base learners. See
#' \link{combinePredictions} for a comprehensive explanation.
#'
#' @slot Mstar Maximum expected loss to be incurred by the \code{baseModels}.
#'
#' @slot wfName workflow key id used to used the right method.
#'
#' @seealso \code{\link{predict}} method for predicting new data using a \code{tseModel}
#' object, \code{\link{forecast}} to forecast new upcoming data points of a time series.
#' \code{\link{SVM}},\code{\link{FFNN}}, \code{\link{RandomForest}} for some examples of
#' base learners implementations.
#'
#' @examples
#' \dontrun{
#' .tse <- tseModel(baseModels,
#'                  preweights,
#'                  target ~.,
#'                  ma.N = NULL,
#'                  embedding.dimension = 30,
#'                  committee.ratio = .2,
#'                  aggregationFUN = "regret",
#'                  Mstar = 30.)
#' }
#' @export
setClass("tseModel",
         slots = c(baseModels = "list",
                   preWeights = "OptionalNumeric",
                   N = "numeric",
                   modelDist = "table",
                   form = "formula",
                   ma.N = "OptionalNumeric",
                   embedding.dimension = "numeric",
                   committee.ratio = "OptionalNumeric",
                   aggregationFUN = "OptionalCharacter",
                   Mstar = "OptionalNumeric",
                   wfName = "character")
)

#' tseModel
#'
#' tseModel is an S4 class that contains the ensemble model. Besides
#' the base learning algorithms -- \code{baseModels} -- tseModel class
#' contains information about other meta-data used to compute predictions
#' for new upcoming data.
#'
#' @param baseModels list comprising the base learners.
#' @param preWeights Normalized relative weights of the base learners according to
#' their performance on the available data.
#' @param form Formula of class formula. Essentially used for reference.
#' @param ma.N The number of periods to average over the Squared Error when
#' computing \strong{MASE} (\emph{Moving Average Squared Error}). This parameter is only
#' used when the function applied to weight the base models is based on
#' MASE. The rationale behind this simple heuristic is that base models are weighted according
#' to their recent performance. Here, recent performance is formalized and quantified as
#' inversely proportional to the models' MASE. The dynamics of the moving average
#' yield a flexibility to the combined model, in the sense that it is self-adaptable when
#' concept drift occurs. Particularly, this parameter \code{ma.N} controls
#' the reactiveness of the learning system to such events.
#' A small value of P leads to greater reactiveness, but also makes the ensemble
#' susceptible to be fooled by outliers. Conversely, higher values of \code{ma.N}
#' lead to greater stability, while sacrificing some responsiveness. This trade-off
#' is known in the literature as the stability-plasticity dilemma (Carpenter et al., 1991).
#' @param embedding.dimension The maximum embedding dimension used to transform the series.
#' @param committee.ratio A numeric value between 0 and 1 representing the ratio
#' of base learners that comprise the committee at each prediction time.
#' If \code{committee.ratio} equals, say, 0.2, at time \emph{t}, the ensemble will the
#' 20\% best base learners up to time \emph{t - 1}.
#' @param aggregationFUN The function name used to combine the base learners. See \link{combinePredictions}
#' for a comprehensive explanation.
#' @param Mstar Maximum expected loss to be incurred by the \code{baseModels}.
#' @seealso \code{\link{predict}} method for predicting new data using a \code{tseModel}
#' object, \code{\link{forecast}} to forecast new upcoming data points of a time series.
#' \code{\link{SVM}},\code{\link{FFNN}}, \code{\link{RandomForest}} for some examples of
#' base learners implementations.
#' @examples
#' \dontrun{
#' .tse <- tseModel(baseModels,
#'                  preweights,
#'                  target ~.,
#'                  ma.N = NULL,
#'                  embedding.dimension = 30,
#'                  committee.ratio = .2,
#'                  aggregationFUN = "regret",
#'                  Mstar = 30.)
#' }
#' @export
tseModel <- function(baseModels,
                     preWeights,
                     form,
                     ma.N = NULL,
                     embedding.dimension,
                     committee.ratio = NULL,
                     aggregationFUN = NULL,
                     Mstar) {
  N <- length(baseModels)
  modelDist <- table(sapply(names(baseModels), function(i) .splitBy_(i)[1]))

  .names <- unique(names(modelDist))
  if ("decisiontree" %in% .names && length(.names) == 1L) {
    wfName <- "SHoW"
  } else if (aggregationFUN %in% c("static-s", "static-w")) {
    wfName <- "SHeW"
  } else {
    wfName <- "DHeW"
  }

  new("tseModel",
      baseModels = baseModels,
      preWeights = preWeights,
      form = form,
      N = N,
      modelDist = modelDist,
      ma.N = ma.N,
      embedding.dimension = embedding.dimension,
      committee.ratio = committee.ratio,
      aggregationFUN = aggregationFUN,
      Mstar = Mstar,
      wfName = wfName)
}
