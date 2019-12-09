setClassUnion("OptionalNumeric", c("numeric", "NULL"))
setClassUnion("OptionalList", c("list","NULL"))


#' committee_set-class
#'
#' Class for a committee set of models. It
#' contains most of the information regarding
#' a \code{\link{constructive_aggregation-class}}
#'
#' @slot committee_set list of subsets
#' @slot form formula
#' @slot specs object of class \code{\link{model_specs-class}}. Contains
#' the parameter setting information for training the
#' base models;
#' @slot out_of_bag Out of bag observations used to compute
#' the subsets
#' @slot recent_series the most recent \code{lambda} observations.
#' @slot lambda window size used to average loss. How far to
#' to back in time.
#' @slot alpha contiguity size. How long should a subset
#' outperform other for it to be considered in the committee
#' @slot aggregate_subsets aggregation approach for the set of subsets.
#' @slot aggregate_hypos final aggregation approach. How should the
#' combined opinions be aggregated.
#'
#' @keywords internal
setClass("committee_set",
         slots = c(committee_set = "list",
                   form = "formula",
                   specs = "model_specs",
                   out_of_bag = "OptionalList",
                   recent_series = "data.frame",
                   lambda = "numeric",
                   alpha = "numeric",
                   aggregate_subsets = "character",
                   aggregate_hypos = "character",
                   covered_regions = "list")
)

#' constructive_aggregation-class
#'
#' Constructive aggregation via out-performance contiguity
#' This method denotes the idea of rearranging a portfolio of
#' models (\strong{base_ensemble}) into different overlapping
#' subsets. These subsets are aggregated (\strong{aggregate_subsets})
#' into combined opinions, forming new models. These models
#' are combined into a final decision through \strong{aggregate_hypos}.
#'
#' @slot base_ensemble object of class \code{\link{base_ensemble-class}}.
#' It contains the base models used that can be used for predicting
#' new data or forecasting future values;
#'
#' @slot committee_set a list of ids -- the individual models
#' of each subset;
#'
#' @slot form formula;
#'
#' @slot specs object of class \code{\link{model_specs-class}}. Contains
#' the parameter setting information for training the
#' base models;
#'
#' @slot lambda window size used to average loss. How far to
#' to back in time.
#'
#' @slot alpha contiguity size. How long should a subset
#' outperform other for it to be considered in the committee
#'
#' @slot recent_series the most recent \code{lambda} observations.
#'
#' @slot out_of_bag Out of bag observations used to compute
#' the subsets
#'
#' @slot aggregate_subsets aggregation approach for the set of subsets.
#'
#' @slot aggregate_hypos final aggregation approach. How should the
#' combined opinions be aggregated.
#'
#' @keywords internal
#'
#' @export
setClass("constructive_aggregation",
         slots = c(base_ensemble = "base_ensemble",
                   committee_set = "OptionalList",
                   form = "formula",
                   specs = "model_specs",
                   lambda = "OptionalNumeric",
                   alpha = "numeric",
                   recent_series = "data.frame",
                   out_of_bag = "OptionalList",
                   aggregate_subsets = "character",
                   aggregate_hypos = "character")
)

#' Build committee set
#'
#' A function for creating the committee set without
#' the final training of the models.
#'
#' @param form formula
#' @param data training data
#' @param specs object of class \code{\link{model_specs-class}}. Contains
#' the parameter setting information for training the
#' base models;
#' @param lambda smoothing window size
#' @param alpha contiguity interval size
#' @param depth depth size how large is the maximum size
#' of the subsets. If NULL, defaults to no. of predictors
#' minus one.
#' @param aggregate_subsets aggregation approach for the
#' set of subsets.
#' @param aggregate_hypos final aggregation approach. How should the
#' combined opinions be aggregated.
#'
#' @export
build_committee_set <-
  function(form,
           data,
           specs,
           lambda,
           alpha,
           depth = NULL,
           aggregate_subsets,
           aggregate_hypos) {

    OOB_data <-
      holdout(
        x = data,
        beta = .7,
        intraining_estimations,
        form = form,
        specs = specs,
        lfun = ae,
        num_cores = 1)

    C <-
      prune_c_outperformance(OOB_data$mloss, lambda, depth)

    C <-
      lapply(C,
             function(x) {
               lapply(x,
                      function(o) sort(o))
             })

    C_star <- prune_c_contiguity(C, alpha)
    covered_regions <- C_star[[2]]
    C_star <- C_star[[1]]

    recent_lambda_k <- recent_lambda_observations(data, lambda)

    methods::new(
      "committee_set",
      committee_set = C_star,
      form = form,
      specs = specs,
      out_of_bag = OOB_data,
      recent_series = recent_lambda_k,
      lambda = lambda,
      alpha = alpha,
      aggregate_subsets = aggregate_subsets,
      aggregate_hypos = aggregate_hypos,
      covered_regions = covered_regions)
  }


#' constructive_aggregation_
#'
#' Simpler way of creating a \code{\link{constructive_aggregation-class}}
#' class object, from \code{\link{committee_set-class}} and
#' \code{\link{base_ensemble-class}}. It assumes that
#' \code{\link{model_specs-class}} are consistent.
#'
#' @param committee_obj \code{\link{committee_set-class}} object
#' @param base_ensemble_obj \code{\link{base_ensemble-class}} object
#'
#'
#' @export
constructive_aggregation_ <-
  function(committee_obj, base_ensemble_obj) {
    methods::new(
      "constructive_aggregation",
      base_ensemble = base_ensemble_obj,
      committee_set = committee_obj@committee_set,
      form = committee_obj@form,
      specs = committee_obj@specs,
      lambda = committee_obj@lambda,
      alpha = committee_obj@alpha,
      recent_series = committee_obj@recent_series,
      out_of_bag = committee_obj@out_of_bag,
      aggregate_subsets = committee_obj@aggregate_subsets,
      aggregate_hypos = committee_obj@aggregate_hypos
    )
  }


#' Constructive aggregation constructor
#'
#' Constructive aggregation via out-performance contiguity
#' This method denotes the idea of rearranging a portfolio of
#' models (\strong{base_ensemble}) into different overlapping
#' subsets. These subsets are aggregated (\strong{aggregate_subsets})
#' into combined opinions, forming new models. These models
#' are combined into a final decision through \strong{aggregate_hypos}.
#'
#' @param form formula
#'
#' @param data training data
#'
#' @param specs object of class \code{\link{model_specs-class}}. Contains
#' the parameter setting information for training the
#' base models;
#' @param lambda smoothing window size
#'
#' @param alpha contiguity interval size
#'
#' @param depth depth size how large is the maximum size
#' of the subsets. If NULL, defaults to no. of predictors
#' minus one.
#'
#' @param aggregate_subsets aggregation approach for the
#' set of subsets.
#' @param aggregate_hypos final aggregation approach. How should the
#' combined opinions be aggregated.
#'
#' @examples
#' specs <- model_specs(
#'   learner = c("bm_svr", "bm_mars"),
#'   learner_pars = list(
#'     bm_glm = list(alpha = c(0, .5, 1)),
#'     bm_svr = list(kernel = c("rbfdot"),
#'                   C = c(1, 3))
#'   )
#' )
#'
#' data("water_consumption")
#' waterc <- embed_timeseries(water_consumption, 5)
#' train <- waterc[1:300, ] # toy size for checks
#' test  <- waterc[301:320, ] # toy size for checks
#'
#' model <- constructive_aggregation(target ~., train, specs, 10,5,NULL,"window_loss","simple")
#' preds <- predict(model, test)
#'
#' @export
constructive_aggregation <-
  function(form,
           data,
           specs,
           lambda = 100,
           alpha = 30,
           depth = NULL,
           aggregate_subsets = "simple",
           aggregate_hypos = "simple") {

    OOB_data <-
      holdout(
        x = data,
        beta = .7,
        intraining_estimations,
        form = form,
        specs = specs,
        lfun = ae,
        num_cores=1)

    C <- prune_c_outperformance(OOB_data$mloss, lambda, depth)

    C <-
      lapply(C,
             function(x) {
               lapply(x,
                      function(o) sort(o))
             })

    C_star <- prune_c_contiguity(C, alpha)

    recent_lambda_k <-
      recent_lambda_observations(data, lambda)

    M <- build_base_ensemble(form, data, specs, 1)

    methods::new(
      "constructive_aggregation",
      base_ensemble = M,
      committee_set = C_star,
      form = form,
      specs = specs,
      lambda = lambda,
      alpha = alpha,
      recent_series = recent_lambda_k,
      out_of_bag = OOB_data,
      aggregate_subsets = aggregate_subsets,
      aggregate_hypos = aggregate_hypos
    )
  }


#' predict method for constructive aggregation
#'
#' Method for predicting new observations using
#' a \code{\link{constructive_aggregation-class}} object
#'
#' @param object a \code{\link{constructive_aggregation-class}} object
#' @param newdata new observations for prediction
#'
#' @keywords internal
#'
#' @export
setMethod("predict",
          signature("constructive_aggregation"),
          function(object, newdata) {
            preds <- hat_info(object, newdata)

            y_hat <-
              combine_committees(preds,
                                 object@aggregate_hypos,
                                 object,
                                 newdata)

            y_hat
          })

setMethod("show",
          signature("constructive_aggregation"),
          function(object) {
            cat("## Constructive Aggregation ##\n\n")

            cat("Base ensemble specs:\n")
            print(object@base_ensemble)
            cat("\n\n")

            cat("Lambda is set to:", object@lambda, "\n")
            cat("Alpha is set to:", object@alpha, "\n")
          })
