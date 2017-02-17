setClassUnion("OptionalNumeric", c("numeric","NULL"))
setClassUnion("OptionalList", c("list","NULL"))
setClassUnion("OptionalDF", c("data.frame","NULL"))

#' Setup base learning models
#'
#' This class sets up the base learning models and respective
#' parameters setting to learn the ensemble.
#'
#' @slot learner Character vector with the base learners
#' to be trained (e.g. \strong{MARS})
#'
#' @slot learner_pars List with parameter setting for the
#' \code{learner}
#'
#' @slot k_size Embedding dimension use to embed the time series
#' into an Euclidean space
#'
#' @slot varying_k An embed split parameter used to encourage
#' data diversity into the ensemble. Each \code{learner} is trained
#' in each setting of \strong{varying_k}. For example, a
#' \strong{varying_k} equal to the \code{k_size} uses the whole
#' embedding dimension; A \strong{varying_k} as \code{c(k_size,
#' k_size / 2)} trains the base models in the whole set of embed
#' cols and in the set of the \strong{varying_k / 2} most
#' recent embeds. Defaults to the whole set (\strong{k_size})
#'
#' @slot varying_w A data split parameter used to encourage
#' data diversity into the ensemble. Similar to \code{varying_k}
#' parameter but related to the embedding vectors instead of the
#' embedding colums. Defaults to the whole data set \code{nrow(data)}
#'
#' @slot data embedded time series of class \code{data.frame}
#'
#' @export
setClass("tseSpecs",
         slots = c(learner = "character",
                   learner_pars = "OptionalList",
                   k_size = "OptionalNumeric",
                   varying_k = "OptionalNumeric",
                   varying_w = "OptionalNumeric",
                   data = "OptionalDF")
)

#' Setup base learning models
#'
#' This class sets up the base learning models and respective
#' parameters setting to learn the ensemble.
#'
#' @param learner Character vector with the base learners
#' to be trained (e.g. \strong{MARS})
#'
#' @param learner_pars List with parameter setting for the
#' \code{learner}
#'
#' @param k_size Embedding dimension use to embed the time series
#' into an Euclidean space
#'
#' @param varying_k An embed split parameter used to encourage
#' data diversity into the ensemble. Each \code{learner} is trained
#' in each setting of \strong{varying_k}. For example, a
#' \strong{varying_k} equal to the \code{k_size} uses the whole
#' embedding dimension; A \strong{varying_k} as \code{c(k_size,
#' k_size / 2)} trains the base models in the whole set of embed
#' cols and in the set of the \strong{varying_k / 2} most
#' recent embeds. Defaults to the whole set (\strong{k_size})
#'
#' @param varying_w A data split parameter used to encourage
#' data diversity into the ensemble. Similar to \code{varying_k}
#' parameter but related to the embedding vectors instead of the
#' embedding colums. Defaults to the whole data set \code{nrow(data)}
#'
#' @param data embedded time series of class \code{data.frame}
#'
#' @export
tseSpecs <- function(learner = "character",
                     learner_pars = "OptionalList",
                     k_size = "OptionalNumeric",
                     varying_k = "OptionalNumeric",
                     varying_w = "OptionalNumeric",
                     data = "OptionalDF") {

  available_learners <- c("baggedtrees", "SVM", "FFNN",
                          "MARS", "Cubist", "RandomForest",
                          "GBM", "GLM", "PPR", "GP", "SAE")

  if (!all(learner %in% available_learners)) {
    stop("All learners must be in list of available learners.\n
         Type ?tseSpecs to check available learners.\n
         User defined learners are not implemented yet.")
  }

  # checks on learner pars list

  if (missing(k_size) || is.null(k_size)) {
    k_size <- get_embedsize(data)
  }

  if (missing(varying_k) || is.null(varying_k)) {
    varying_k <- k_size
  }

  if (missing(varying_w) || is.null(varying_w)) {
    varying_w <- nrow(data)
  }

  new("tseSpecs",
      learner = learner,
      learner_pars = learner_pars,
      k_size = k_size,
      varying_k = varying_k,
      varying_w = varying_w,
      data = data)
  }


setMethod("show",
          signature("tseSpecs"),
          function(object) {
            clearner <- paste(object@learner, collapse = ", ")

            cat("Ensemble setup with the individual models:", clearner,".\n\n")

            clpars <- paste(names(object@learner_pars), collapse = ", ")

            cat("With parameters set up for learners: ", clpars, "\n")
            cat("Others will be set up by default\n\n")
            cat("Embedding size set up to:", object@k_size, "\n")
            cat("Embed window splitting set to:", object@varying_k, "\n")
            cat("Data window splitting set to:", object@varying_w, "\n")
          }
)
