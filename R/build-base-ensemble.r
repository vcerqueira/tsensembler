#' Wrapper for creating an ensemble
#'
#' Using the parameter specifications from
#' \code{\link{model_specs-class}}, this function trains
#' a set of regression models.
#'
#' @param form formula;
#' @param data data.frame for training the predictive models;
#' @param specs object of class \code{\link{model_specs-class}}. Contains the information
#' about the parameter setting of the models to train.
#' @param num_cores number of cores
#'
#'
#' @return An S4 class with the following slots:
#' \strong{base_models}, a list containing the trained models;
#' \strong{pre_weights}, a numeric vector describing the weights
#' of the base models according to their performance in the training
#' data; and \strong{colnames}, the column names of the data, used for
#' reference.
#'
#' @examples
#' data("water_consumption")
#' dataset <- embed_timeseries(water_consumption, 5)
#' specs <- model_specs(c("bm_ppr","bm_svr"), NULL)
#' M <- build_base_ensemble(target ~., dataset, specs, 1)
#'
#' @export
build_base_ensemble <-
  function(form, data, specs, num_cores=1) {
    M <- learning_base_models(data, form, specs, num_cores=num_cores)

    base_ensemble(base_models = M$base_model,
                  pre_weights =  M$preweights,
                  form = form,
                  colnames = colnames(data))
  }

#' Training the base models of an ensemble
#'
#' This function uses \emph{train} to build a set
#' of predictive models, according to \emph{specs}
#'
#' @param train training set to build the predictive models;
#' @param form formula;
#' @param specs object of class \code{\link{model_specs-class}}
#' @param num_cores A numeric value to specify the number of cores used to
#' train base and meta models. num_cores = 1
#' leads to sequential training of models. num_cores > 1
#' splits the training of the base models across num_cores cores.
#'
#' @seealso \code{\link{build_base_ensemble}}.
#'
#' @examples
#' data("water_consumption")
#' dataset <- embed_timeseries(water_consumption, 5)
#' specs <- model_specs(c("bm_ppr","bm_svr"), NULL)
#' M <- build_base_ensemble(target ~., dataset, specs, 1)
#'
#' @return A series of predictive models (\code{base_model}), and
#' the weights of the models computed in the training
#' data (\code{preweights}).
learning_base_models <-
  function(train, form, specs, num_cores) {
    if (length(num_cores) > 1 || !is.numeric(num_cores)) {
      stop("Please specify a numeric value for num_cores. num_cores = 1
           leads to sequential training of models. num_cores > 1
           splits the training of the base models across num_cores cores.")
    }

    if (is.null(num_cores)) num_cores <- 1

    Y_tr <- get_y(train, form)

    learner <- specs@learner
    lpars <- specs@learner_pars

    cat("\nTraining the base models...\n")
    if (num_cores > 1) {
      cl <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl)
      `%partrain%` <- `%dopar%`
    } else {
      `%partrain%` <- `%do%`
    }

    o <- NULL
    base_model <-
      foreach::foreach(o = seq_along(learner),
                       .packages = "tsensembler") %partrain% {
                         cat("\n\n",learner[o],"\n")

                         do.call(learner[o], c(list(form, train, lpars)))
                       }

    pre_weights <-
      lapply(base_model,
           function(o) {
             compute_predictions(o, form, train)
           })

    base_model <- unlist(base_model, recursive = FALSE)
    pre_weights  <- unlist(pre_weights, recursive = FALSE)

    W <- vapply(pre_weights, function(o) {
      rmse(Y_tr, o)
    }, numeric(1), USE.NAMES = FALSE)
    W <- model_weighting(W, trans = "linear")

    list(base_model = base_model, preweights = W)
  }

setMethod("show",
          signature("base_ensemble"),
          function(object) {
            cat("Time Series Ensemble Model\n")
            cat("Ensemble composed by", object@N, "base Models.\n")
            cat("The distribution of the base Models is the following:\n")
            print(object@model_distribution)
            cat("The formula is:", deparse(object@form))
          })
