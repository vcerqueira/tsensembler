setClassUnion("OptionalNumeric", c("numeric", "NULL"))
setClassUnion("OptionalList", c("list","NULL"))
setClassUnion("OptionalCharacter", c("character","NULL"))

#' base_ensemble-class
#'
#' \strong{base_ensemble} is a S4 class that contains the base models
#' comprising the ensemble. Besides the base learning algorithms --
#' \code{base_models} -- base_ensemble class contains information
#' about other meta-data used to compute predictions for new upcoming data.
#'
#' @slot base_models a list comprising the base models;
#'
#' @slot pre_weights Normalized relative weights of the base learners according to
#' their performance on the available data;
#'
#' @slot form formula;
#'
#' @slot colnames names of the columns of the data used to train the \strong{base_models};
#'
#' @slot N number of base models;
#'
#' @slot model_distribution base learner distribution with respect to the type of learner.
#' That is, the number of Decision Trees, SVMs, etc.
#'
#' @keywords internal
#'
#' @export
setClass("base_ensemble",
         slots = c(base_models = "list",
                   pre_weights = "OptionalNumeric",
                   form = "formula",
                   colnames = "OptionalCharacter",
                   N = "numeric",
                   model_distribution = "numeric")
)

#' base_ensemble
#'
#' \strong{base_ensemble} is a S4 class that contains the base models
#' comprising the ensemble. Besides the base learning algorithms --
#' \code{base_models} -- base_ensemble class contains information
#' about other meta-data used to compute predictions for new upcoming data.
#'
#' @param base_models a list comprising the base models;
#' @param pre_weights normalized relative weights of the base learners according to
#' their performance on the available data;
#' @param form formula;
#' @param colnames names of the columns of the data used to train the \strong{base_models};
#'
#' @export
base_ensemble <-
  function(base_models, pre_weights, form, colnames) {
    N <- length(base_models)
    mnames <- sapply(names(base_models),
                     function(i) split_by_(i)[1])


    model_distribution <- table(mnames)
    unames <- names(model_distribution)

    model_distribution <- as.vector(model_distribution)
    names(model_distribution) <- unames


    methods::new(
      "base_ensemble",
      base_models = base_models,
      pre_weights = pre_weights,
      form = form,
      colnames = colnames,
      N = N,
      model_distribution = model_distribution
    )
  }


#' Setup base learning models
#'
#' This class sets up the base learning models and respective
#' parameters setting to learn the ensemble.
#'
#' @slot learner character vector with the base learners
#' to be trained. Currently available models are:
#' \describe{
#' \item{\strong{bm_gaussianprocess}}{Gaussian Process models, from the
#' \strong{kernlab} package. See \code{\link[kernlab]{gausspr}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_gaussianprocess}} for the function implementation.}
#'
#' \item{\strong{bm_ppr}}{Projection Pursuit Regression models, from the
#' \strong{stats} package. See \code{\link[stats]{ppr}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_ppr}} for the function implementation.}
#'
#' \item{\strong{bm_glm}}{Generalized Linear Models, from the
#' \strong{glmnet} package. See \code{\link[glmnet]{glmnet}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_glm}} for the function implementation.}
#'
#' \item{\strong{bm_gbm}}{Generalized Boosted Regression models, from the
#' \strong{gbm} package. See \code{\link[gbm]{gbm}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_gbm}} for the function implementation.}
#'
#' \item{\strong{bm_randomforest}}{Random Forest models, from the
#' \strong{ranger} package. See \code{\link[ranger]{ranger}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_randomforest}} for the function implementation.}
#'
#' \item{\strong{bm_cubist}}{M5 tree models, from the
#' \strong{Cubist} package. See \code{\link[Cubist]{cubist}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_cubist}} for the function implementation.}
#'
#' \item{\strong{bm_mars}}{Multivariate Adaptive Regression Splines models, from the
#' \strong{earth} package. See \code{\link[earth]{earth}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_mars}} for the function implementation.}
#'
#' \item{\strong{bm_svr}}{Support Vector Regression models, from the
#' \strong{kernlab} package. See \code{\link[kernlab]{ksvm}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_svr}} for the function implementation.}
#'
#' \item{\strong{bm_ffnn}}{Feedforward Neural Network models, from the
#' \strong{nnet} package. See \code{\link[nnet]{nnet}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_ffnn}} for the function implementation.}
#'
#' \item{\strong{bm_pls_pcr}}{Partial Least Regression and Principal
#' Component Regression models, from the \strong{pls} package.
#' See \code{\link[pls]{mvr}} for a complete description
#' and possible parametrization. See \code{\link{bm_pls_pcr}}
#' for the function implementation.}
#' }
#'
#' @slot learner_pars a list with parameter setting for the
#' \strong{learner}. For each model, a inner list should be created
#' with the specified parameters.
#'
#' Check each implementation to see the possible
#' variations of parameters (also examplified below).
#'
#' @examples
#' # A PPR model and a GLM model with default parameters
#' model_specs(learner = c("bm_ppr", "bm_glm"), learner_pars = NULL)
#'
#'
#' # A PPR model and a SVR model. The listed parameters are combined
#' # with a cartesian product.
#' # With these specifications an ensemble with 6 predictive base
#' # models will be created. Two PPR models, one with 2 nterms
#' # and another with 4; and 4 SVR models, combining the kernel
#' # and C parameters.
#' specs <- model_specs(
#'  c("bm_ppr", "bm_svr"),
#'  list(bm_ppr = list(nterms = c(2, 4)),
#'       bm_svr = list(kernel = c("vanilladot", "polydot"), C = c(1,5)))
#' )
#'
#' # All parameters currently available (parameter values can differ)
#' model_specs(
#'  learner = c("bm_ppr", "bm_svr", "bm_randomforest",
#'              "bm_gaussianprocess", "bm_cubist", "bm_glm",
#'              "bm_gbm", "bm_pls_pcr", "bm_ffnn", "bm_mars"
#'          ),
#'  learner_pars = list(
#'     bm_ppr = list(
#'        nterms = c(2,4),
#'        sm.method = "supsmu"
#'      ),
#'     bm_svr = list(
#'        kernel = "rbfdot",
#'        C = c(1,5),
#'        epsilon = .01
#'      ),
#'     bm_glm = list(
#'        alpha = c(1, 0)
#'      ),
#'     bm_randomforest = list(
#'        num.trees = 500
#'      ),
#'     bm_gbm = list(
#'        interaction.depth = 1,
#'        shrinkage = c(.01, .005),
#'        n.trees = c(100)
#'      ),
#'     bm_mars = list(
#'        nk = 15,
#'        degree = 3,
#'        thresh = .001
#'      ),
#'     bm_ffnn = list(
#'        size = 30,
#'        decay = .01
#'      ),
#'     bm_pls_pcr = list(
#'        method = c("kernelpls", "simpls", "cppls")
#'      ),
#'     bm_gaussianprocess = list(
#'        kernel = "vanilladot",
#'        tol = .01
#'      ),
#'     bm_cubist = list(
#'        committees = 50,
#'        neighbors = 0
#'      )
#'   )
#' )
#'
#' @export
setClass("model_specs",
         slots = c(learner = "character", learner_pars = "OptionalList")
)

#' Setup base learning models
#'
#' This class sets up the base learning models and respective
#' parameters setting to learn the ensemble.
#'
#' @param learner character vector with the base learners
#' to be trained. Currently available models are:
#' \describe{
#' \item{\strong{bm_gaussianprocess}}{Gaussian Process models, from the
#' \strong{kernlab} package. See \code{\link[kernlab]{gausspr}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_gaussianprocess}} for the function implementation.}
#'
#' \item{\strong{bm_ppr}}{Projection Pursuit Regression models, from the
#' \strong{stats} package. See \code{\link[stats]{ppr}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_ppr}} for the function implementation.}
#'
#' \item{\strong{bm_glm}}{Generalized Linear Models, from the
#' \strong{glmnet} package. See \code{\link[glmnet]{glmnet}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_glm}} for the function implementation.}
#'
#' \item{\strong{bm_gbm}}{Generalized Boosted Regression models, from the
#' \strong{gbm} package. See \code{\link[gbm]{gbm}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_gbm}} for the function implementation.}
#'
#' \item{\strong{bm_randomforest}}{Random Forest models, from the
#' \strong{ranger} package. See \code{\link[ranger]{ranger}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_randomforest}} for the function implementation.}
#'
#' \item{\strong{bm_cubist}}{M5 tree models, from the
#' \strong{Cubist} package. See \code{\link[Cubist]{cubist}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_cubist}} for the function implementation.}
#'
#' \item{\strong{bm_mars}}{Multivariate Adaptive Regression Splines models, from the
#' \strong{earth} package. See \code{\link[earth]{earth}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_mars}} for the function implementation.}
#'
#' \item{\strong{bm_svr}}{Support Vector Regression models, from the
#' \strong{kernlab} package. See \code{\link[kernlab]{ksvm}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_svr}} for the function implementation.}
#'
#' \item{\strong{bm_ffnn}}{Feedforward Neural Network models, from the
#' \strong{nnet} package. See \code{\link[nnet]{nnet}}
#' for a complete description and possible parametrization. See
#' \code{\link{bm_ffnn}} for the function implementation.}
#'
#' \item{\strong{bm_pls_pcr}}{Partial Least Regression and Principal
#' Component Regression models, from the \strong{pls} package.
#' See \code{\link[pls]{mvr}} for a complete description
#' and possible parametrization. See \code{\link{bm_pls_pcr}}
#' for the function implementation.}
#' }
#'
#' @param learner_pars a list with parameter setting for the
#' \strong{learner}. For each model, a inner list should be created
#' with the specified parameters.
#'
#' Check each implementation to see the possible
#' variations of parameters (also examplified below).
#'
#' @examples
#' # A PPR model and a GLM model with default parameters
#' model_specs(learner = c("bm_ppr", "bm_glm"), learner_pars = NULL)
#'
#'
#' # A PPR model and a SVR model. The listed parameters are combined
#' # with a cartesian product.
#' # With these specifications an ensemble with 6 predictive base
#' # models will be created. Two PPR models, one with 2 nterms
#' # and another with 4; and 4 SVR models, combining the kernel
#' # and C parameters.
#' specs <- model_specs(
#'  c("bm_ppr", "bm_svr"),
#'  list(bm_ppr = list(nterms = c(2, 4)),
#'       bm_svr = list(kernel = c("vanilladot", "polydot"), C = c(1,5)))
#' )
#'
#' # All parameters currently available (parameter values can differ)
#' model_specs(
#'  learner = c("bm_ppr", "bm_svr", "bm_randomforest",
#'              "bm_gaussianprocess", "bm_cubist", "bm_glm",
#'              "bm_gbm", "bm_pls_pcr", "bm_ffnn", "bm_mars"
#'          ),
#'  learner_pars = list(
#'     bm_ppr = list(
#'        nterms = c(2,4),
#'        sm.method = "supsmu"
#'      ),
#'     bm_svr = list(
#'        kernel = "rbfdot",
#'        C = c(1,5),
#'        epsilon = .01
#'      ),
#'     bm_glm = list(
#'        alpha = c(1, 0)
#'      ),
#'     bm_randomforest = list(
#'        num.trees = 500
#'      ),
#'     bm_gbm = list(
#'        interaction.depth = 1,
#'        shrinkage = c(.01, .005),
#'        n.trees = c(100)
#'      ),
#'     bm_mars = list(
#'        nk = 15,
#'        degree = 3,
#'        thresh = .001
#'      ),
#'     bm_ffnn = list(
#'        size = 30,
#'        decay = .01
#'      ),
#'     bm_pls_pcr = list(
#'        method = c("kernelpls", "simpls", "cppls")
#'      ),
#'     bm_gaussianprocess = list(
#'        kernel = "vanilladot",
#'        tol = .01
#'      ),
#'     bm_cubist = list(
#'        committees = 50,
#'        neighbors = 0
#'      )
#'   )
#' )
#'
#' @export
model_specs <-
  function(learner, learner_pars = NULL) {
    .available_models <-
      c("bm_gaussianprocess",
        "bm_ppr",
        "bm_glm",
        "bm_gbm",
        "bm_randomforest",
        "bm_ffnn",
        "bm_svr",
        "bm_mars",
        "bm_cubist",
        "bm_pls_pcr"
        )

    if (!all(learner %in% .available_models))
      stop("One or more invalid base models.", call. = FALSE)

    if (!is.null(learner_pars)) {
      if (!all(names(learner_pars) %in% .available_models)) {
        warning("Some model parameter name badly specified.
                The parameter list must have the same names as the models.
                Check ?model_specs for the available models and some examples.")
      }

      names_pars <- names(learner_pars)
      for (model in learner)
        if (is_model_in_pars(model,
                             learner,
                             names_pars)) {
          are_pars_valid(model, learner_pars)
        }
    }

    methods::new("model_specs",
                 learner = learner,
                 learner_pars = learner_pars)
  }

setMethod("show",
          signature("model_specs"),
          function(object) {
            mnames <- paste(object@learner, collapse = ",\n ")

            cat("Ensemble setup with the individual models:\n", mnames,"\n\n")

            if (is.null(object@learner_pars))
              cat("Parameter setup will be set with default values.\n")
            else {
              m_tbl <- lapply(object@learner_pars, expand.grid)

              null.pars <- !(object@learner %in% names(m_tbl))

              set_up_models <- object@learner[!null.pars]
              par_models <- paste(set_up_models, collapse = ", ")
              cat(par_models,
                  "models will have the following parameters set up:\n\n")
              for (m in set_up_models) {
                cat(m, ": ")
                cat(paste(names(object@learner_pars[[m]]), collapse = ", "), "\n")
              }

              null_par_models <-
                paste(object@learner[null.pars], collapse = ", ")
              cat(null_par_models,
                  "models will be set up with default parameters.\n")

              N <- sum(sapply(m_tbl, nrow)) + sum(null.pars)

              cat("\nWith these specs, the ensemble size is:", N)
            }
          }
)


#' Arbitrated Dynamic Ensemble
#'
#' Arbitrated Dynamic Ensemble (ADE) is an ensemble approach
#' for adaptively combining forecasting models. A metalearning
#' strategy is used that specializes base models
#' across the time series. Each meta-learner is specifically
#' designed to model how apt its base counterpart is to make
#' a prediction for a given test example. This is accomplished
#' by analysing how the error incurred by a given learning model
#' relates to the characteristics of the data. At test time,
#' the base-learners are weighted according to their degree
#' of competence in the input observation, estimated by the
#' predictions of the meta-learners.
#'
#' @slot base_ensemble object of class \code{\link{base_ensemble-class}}.
#' It contains the base models used that can be used for predicting
#' new data or forecasting future values;
#'
#' @slot meta_model a list containing the meta models, one for
#' each base model. The meta-models are random forests;
#'
#' @slot form formula;
#'
#' @slot specs object of class \code{\link{model_specs-class}}. Contains
#' the parameter setting information for training the
#' base models;
#'
#' @slot lambda window size. Number of observations to compute
#' the recent performance of the base models, according to the
#' committee ratio \strong{omega}. Essentially, the top \emph{omega}
#' models are selected and weighted at each prediction instance, according
#' to their performance in the last \emph{lambda} observations.
#' Defaults to 50 according to empirical experiments;
#'
#' @slot omega committee ratio size. Essentially, the top \emph{omega} * 100
#' percent of models are selected and weighted at each prediction instance, according
#' to their performance in the last \emph{lambda} observations.
#' Defaults to .5 according to empirical experiments;
#'
#' @slot select_best Logical. If true, at each prediction time,
#' a single base model is picked to make a prediction. The picked
#' model is the one that has the lowest loss prediction from
#' the meta models. Defaults to FALSE;
#'
#' @slot recent_series the most recent \code{lambda} observations.
#'
#' @references Cerqueira, Vitor; Torgo, Luis; Pinto, Fabio;
#' and Soares, Carlos. "Arbitrated Ensemble for Time Series
#' Forecasting" to appear at: Joint European Conference on Machine Learning and
#' Knowledge Discovery in Databases. Springer International
#' Publishing, 2017.
#'
#' V. Cerqueira, L. Torgo, and C. Soares, “Arbitrated ensemble for
#' solar radiation forecasting,” in International Work-Conference on
#' Artificial Neural Networks. Springer, Cham, 2017, pp. 720–732
#'
#' @seealso \code{\link{model_specs-class}} for setting up the ensemble parameters
#' for an \strong{ADE} model; \code{\link{forecast}} for the forecasting method
#' that uses an \strong{ADE} model for forecasting future values;
#' \code{\link{predict}} for the method that predicts new held out observations;
#' \code{\link{update_weights}} for the method used to update the
#' weights of an \strong{ADE} model between successive predict or forecast calls;
#' \code{\link{update_ade_meta}} for updating (retraining) the meta models
#' of an \strong{ADE} model; \code{\link{update_base_models}} for
#' the updating (retraining) the base models of an \strong{ADE} ensemble (and respective
#' weights); \code{\link{ade_hat-class}} for the object that results from
#' predicting with an \strong{ADE} model; and \code{\link{update_ade}} to update an ADE
#' model, combining functions \strong{update_base_models}, \strong{update_meta_ade}, and
#' \strong{update_weights}.
#'
#' @examples
#' \dontrun{
#' specs <- model_specs(
#'   learner = c("bm_ppr", "bm_glm", "bm_svr", "bm_mars"),
#'   learner_pars = list(
#'     bm_glm = list(alpha = c(0, .5, 1)),
#'     bm_svr = list(kernel = c("rbfdot", "polydot"),
#'                   C = c(1, 3)),
#'     bm_ppr = list(nterms = 4)
#'   )
#' )
#'
#' data("water_consumption")
#' train <- embed_timeseries(water_consumption, 5)
#'
#' model <- ADE(target ~., train, specs)
#' }
#'
#' @export
setClass("ADE",
         slots = c(base_ensemble = "base_ensemble",
                   meta_model = "list",
                   form = "formula",
                   specs = "model_specs",
                   lambda = "numeric",
                   omega = "OptionalNumeric",
                   select_best = "logical",
                   recent_series = "data.frame")
)

#' Arbitrated Dynamic Ensemble
#'
#' Arbitrated Dynamic Ensemble (ADE) is an ensemble approach
#' for adaptively combining forecasting models. A metalearning
#' strategy is used that specializes base models
#' across the time series. Each meta-learner is specifically
#' designed to model how apt its base counterpart is to make
#' a prediction for a given test example. This is accomplished
#' by analysing how the error incurred by a given learning model
#' relates to the characteristics of the data. At test time,
#' the base-learners are weighted according to their degree
#' of competence in the input observation, estimated by the
#' predictions of the meta-learners.
#'
#' @param form formula;
#'
#' @param data data to train the base models
#'
#' @param specs object of class \code{\link{model_specs-class}}. Contains
#' the parameter setting information for training the
#' base models;
#'
#' @param lambda window size. Number of observations to compute
#' the recent performance of the base models, according to the
#' committee ratio \strong{omega}. Essentially, the top \emph{omega}
#' models are selected and weighted at each prediction instance, according
#' to their performance in the last \emph{lambda} observations.
#' Defaults to 50 according to empirical experiments;
#'
#' @param omega committee ratio size. Essentially, the top \emph{omega} * 100
#' percent of models are selected and weighted at each prediction instance, according
#' to their performance in the last \emph{lambda} observations.
#' Defaults to .5 according to empirical experiments;
#'
#' @param select_best Logical. If true, at each prediction time,
#' a single base model is picked to make a prediction. The picked
#' model is the one that has the lowest loss prediction from
#' the meta models. Defaults to FALSE;
#'
#' @references Cerqueira, Vitor; Torgo, Luis; Pinto, Fabio;
#' and Soares, Carlos. "Arbitrated Ensemble for Time Series
#' Forecasting" to appear at: Joint European Conference on Machine Learning and
#' Knowledge Discovery in Databases. Springer International
#' Publishing, 2017.
#'
#' V. Cerqueira, L. Torgo, and C. Soares, “Arbitrated ensemble for
#' solar radiation forecasting,” in International Work-Conference on
#' Artificial Neural Networks. Springer, Cham, 2017, pp. 720–732
#'
#' @seealso \code{\link{model_specs-class}} for setting up the ensemble parameters
#' for an \strong{ADE} model; \code{\link{forecast}} for the forecasting method
#' that uses an \strong{ADE} model for forecasting future values;
#' \code{\link{predict}} for the method that predicts new held out observations;
#' \code{\link{update_weights}} for the method used to update the
#' weights of an \strong{ADE} model between successive predict or forecast calls;
#' \code{\link{update_ade_meta}} for updating (retraining) the meta models
#' of an \strong{ADE} model; \code{\link{update_base_models}} for
#' the updating (retraining) the base models of an \strong{ADE} ensemble (and respective
#' weights); \code{\link{ade_hat-class}} for the object that results from
#' predicting with an \strong{ADE} model; and \code{\link{update_ade}} to update an ADE
#' model, combining functions \strong{update_base_models}, \strong{update_meta_ade}, and
#' \strong{update_weights}.
#'
#' @examples
#' \dontrun{
#' specs <- model_specs(
#'  learner = c("bm_svr", "bm_glm", "bm_mars"),
#'  learner_pars = NULL
#' )
#'
#' data("water_consumption")
#' train <- embed_timeseries(water_consumption, 5)
#'
#' model <- ADE(target ~., train, specs)
#' }
#'
#' @export
"ADE" <-
  function(form,
           data,
           specs,
           lambda = 50,
           omega = .5,
           select_best = FALSE) {

    if (select_best && is.numeric(omega))
      warning(
        "ADE setup with both selection of best learner and
        a committee (\"omega\" parameter).
        \"omega\" parameter will be ignored."
        )

    if (!is.null(omega))
      if (omega >= 1 | omega <= 0)
        stop("\"omega\" parameter must be a double between 0 and 1",
             call. = FALSE)

    if (lambda < 1 | lambda > nrow(data))
      stop("\"lambda\" parameter must be a positive integer < nrow(data)",
           call. = FALSE)

    M <-
      train_ade(
        form = form,
        train = data,
        specs = specs,
        lambda = lambda
      )

    methods::new(
      "ADE",
      base_ensemble = M$base_ensemble,
      meta_model = M$meta_model,
      form = form,
      specs = specs,
      lambda = lambda,
      omega = omega,
      select_best = select_best,
      recent_series = M$recent_series
    )
  }


#' Predictions by an ADE ensemble
#'
#' Predictions produced by a \code{\link{ADE-class}} object.
#' It contains \strong{y_hat}, the combined predictions,
#' \strong{Y_hat}, the predictions of each base model,
#' \strong{Y_committee}, the base models used for prediction
#' at each time point, and \strong{E_hat}, the loss predictions
#' by each meta-model.
#'
#' @slot y_hat combined predictions of the ensemble
#' \code{\link{ADE-class}}. A numeric vector;
#'
#' @slot Y_hat a matrix containing the predictions made by
#' individual models;
#'
#' @slot Y_committee a list describing the models selected for
#' predictions at each time point (according to \strong{lambda}
#' and \strong{omega});
#'
#' @slot E_hat predictions of error of each base model, estimated
#' by their respective meta model associate;
#'
#' @seealso \code{\link{ADE}} for generating an ADE ensemble.
#'
#' @family ensemble predictions
#'
#' @export
setClass("ade_hat",
         slots = c(y_hat = "numeric",
                   Y_hat = "data.frame",
                   Y_committee = "OptionalList",
                   E_hat = "data.frame")
)



#' Predictions by an ADE ensemble
#'
#' Predictions produced by a \code{\link{ADE-class}} object.
#' It contains \strong{y_hat}, the combined predictions,
#' \strong{Y_hat}, the predictions of each base model,
#' \strong{Y_committee}, the base models used for prediction
#' at each time point, and \strong{E_hat}, the loss predictions
#' by each meta-model.
#'
#' @param y_hat combined predictions of the ensemble
#' \code{\link{ADE}}. A numeric vector;
#'
#' @param Y_hat a matrix containing the predictions made by
#' individual models;
#'
#' @param Y_committee a list describing the models selected for
#' predictions at each time point (according to \strong{lambda}
#' and \strong{omega});
#'
#' @param E_hat predictions of error of each base model, estimated
#' by their respective meta model associate;
#'
#' @family ensemble predictions
#'
#' @seealso \code{\link{ADE-class}} for generating an ADE ensemble.
#'
#' @export
ade_hat <- function(y_hat, Y_hat, Y_committee, E_hat) {

  if (missing(Y_committee)) Y_committee <- NULL

  methods::new(
    "ade_hat",
    y_hat = y_hat,
    Y_hat = Y_hat,
    Y_committee = Y_committee,
    E_hat = E_hat
  )
}


#' Dynamic Ensemble for Time Series
#'
#' A Dynamic Ensemble for Time Series (DETS). The DETS ensemble
#' method we present settles on individually pre-trained models
#' which are dynamically combined at run-time to make a prediction.
#' The combination rule is reactive to changes in the environment,
#' rendering an online combined model. The main properties of the ensemble
#' are:
#' \describe{
#'  \item{heterogeneity}{Heterogeneous ensembles are those
#'  comprised of different types of base learners. By employing
#'  models that follow different learning strategies, use different
#'  features and/or data observations we expect that individual
#'  learners will disagree with each other, introducing a natural
#'  diversity into the ensemble that helps in handling different
#'  dynamic regimes in a time series forecasting setting;}
#'  \item{responsiveness}{We promote greater responsiveness of
#'  heterogeneous ensembles in time series tasks by making the
#'  aggregation of their members' predictions time-dependent.
#'  By tracking the loss of each learner over time, we weigh
#'  the predictions of individual learners according to their
#'  recent performance using a non-linear function. This strategy
#'  may be advantageous for better detecting regime changes and
#'  also to quickly adapt the ensemble to new regimes.}
#' }
#'
#' @slot base_ensemble object of class \code{\link{base_ensemble-class}}.
#' It contains the base models used that can be used for predicting
#' new data or forecasting future values;
#'
#' @slot form formula;
#'
#' @slot specs object of class \code{\link{model_specs-class}}. Contains
#' the parameter setting information for training the
#' base models;
#'
#' @slot lambda window size. Number of observations to compute
#' the recent performance of the base models, according to the
#' committee ratio \strong{omega}. Essentially, the top \emph{omega}
#' models are selected and weighted at each prediction instance, according
#' to their performance in the last \emph{lambda} observations.
#' Defaults to 50 according to empirical experiments;
#'
#' @slot omega committee ratio size. Essentially, the top \emph{omega}
#' models are selected and weighted at each prediction instance, according
#' to their performance in the last \emph{lambda} observations.
#' Defaults to .5 according to empirical experiments;
#'
#' @slot recent_series the most recent \code{lambda} observations.
#'
#' @references Cerqueira, Vitor; Torgo, Luis; Oliveira, Mariana,
#' and Bernhard Pfahringer. "Dynamic and Heterogeneous Ensembles
#' for Time Series Forecasting." Data Science and Advanced
#' Analytics (DSAA), 2017 IEEE International Conference on. IEEE, 2017.
#'
#' @seealso \code{\link{model_specs-class}} for setting up the ensemble parameters
#' for an \strong{DETS} model; \code{\link{forecast}} for the method
#' that uses an \strong{DETS} model for forecasting future values;
#' \code{\link{predict}} for the method that predicts new held out observations;
#' \code{\link{update_weights}} for the method used to update the
#' weights of an \strong{DETS} model between successive predict or forecast calls;
#' \code{\link{update_base_models}} for the updating (retraining)
#' the base models of an \strong{DETS} ensemble (and respective
#' weights); and \code{\link{dets_hat-class}} for the object that results from
#' predicting with an \strong{DETS} model.
#'
#' @examples
#' specs <- model_specs(
#'  learner = c("bm_randomforest", "bm_ppr", "bm_mars"),
#'  learner_pars = NULL
#' )
#'
#' data("water_consumption")
#' train <- embed_timeseries(water_consumption, 5)
#'
#' model <- DETS(target ~., train, specs, lambda = 30, omega = .2)
#'
#' @export
setClass("DETS",
         slots = c(base_ensemble = "base_ensemble",
                   form = "formula",
                   specs = "model_specs",
                   lambda = "OptionalNumeric",
                   omega = "OptionalNumeric",
                   recent_series = "data.frame")
)

#' Dynamic Ensemble for Time Series
#'
#' A Dynamic Ensemble for Time Series (DETS). The DETS ensemble
#' method we present settles on individually pre-trained models
#' which are dynamically combined at run-time to make a prediction.
#' The combination rule is reactive to changes in the environment,
#' rendering an online combined model. The main properties of the ensemble
#' are:
#' \describe{
#'  \item{heterogeneity}{Heterogeneous ensembles are those
#'  comprised of different types of base learners. By employing
#'  models that follow different learning strategies, use different
#'  features and/or data observations we expect that individual
#'  learners will disagree with each other, introducing a natural
#'  diversity into the ensemble that helps in handling different
#'  dynamic regimes in a time series forecasting setting;}
#'  \item{responsiveness}{We promote greater responsiveness of
#'  heterogeneous ensembles in time series tasks by making the
#'  aggregation of their members' predictions time-dependent.
#'  By tracking the loss of each learner over time, we weigh
#'  the predictions of individual learners according to their
#'  recent performance using a non-linear function. This strategy
#'  may be advantageous for better detecting regime changes and
#'  also to quickly adapt the ensemble to new regimes.}
#' }
#'
#' @param form formula;
#'
#' @param data data frame to train the base models;
#'
#' @param specs object of class \code{\link{model_specs-class}}. Contains
#' the parameter setting information for training the
#' base models;
#'
#' @param lambda window size. Number of observations to compute
#' the recent performance of the base models, according to the
#' committee ratio \strong{omega}. Essentially, the top \emph{omega}
#' models are selected and weighted at each prediction instance, according
#' to their performance in the last \emph{lambda} observations.
#' Defaults to 50 according to empirical experiments;
#'
#' @param omega committee ratio size. Essentially, the top \emph{omega}
#' models are selected and weighted at each prediction instance, according
#' to their performance in the last \emph{lambda} observations.
#' Defaults to .5 according to empirical experiments;
#'
#' @references Cerqueira, Vitor; Torgo, Luis; Oliveira, Mariana,
#' and Bernhard Pfahringer. "Dynamic and Heterogeneous Ensembles
#' for Time Series Forecasting." Data Science and Advanced
#' Analytics (DSAA), 2017 IEEE International Conference on. IEEE, 2017.
#'
#' @seealso \code{\link{model_specs-class}} for setting up the ensemble parameters
#' for an \strong{DETS} model; \code{\link{forecast}} for the method
#' that uses an \strong{DETS} model for forecasting future values;
#' \code{\link{predict}} for the method that predicts new held out observations;
#' \code{\link{update_weights}} for the method used to update the
#' weights of an \strong{DETS} model between successive predict or forecast calls;
#' \code{\link{update_base_models}} for the updating (retraining)
#' the base models of an \strong{DETS} ensemble (and respective
#' weights); and \code{\link{dets_hat-class}} for the object that results from
#' predicting with an \strong{DETS} model.
#'
#' @examples
#' specs <- model_specs(
#'  learner = c("bm_randomforest", "bm_ppr", "bm_mars"),
#'  learner_pars = NULL
#' );
#'
#' data("water_consumption");
#' train <- embed_timeseries(water_consumption, 5);
#'
#' model <- DETS(target ~., train, specs, lambda = 30, omega = .2)
#'
#' @export
"DETS" <-
  function(form,
           data,
           specs,
           lambda = 50,
           omega = .5) {

    M <- build_base_ensemble(form, data, specs)

    recent_lambda_k <- recent_lambda_observations(data, lambda)

    methods::new(
      "DETS",
      base_ensemble = M,
      form = form,
      specs = specs,
      lambda = lambda,
      omega = omega,
      recent_series = recent_lambda_k
    )
  }

#' Predictions by an DETS ensemble
#'
#' @slot y_hat combined predictions of the ensemble
#' \code{\link{DETS-class}}. A numeric vector;
#'
#' @slot Y_hat a matrix containing the predictions made by
#' individual models;
#'
#' @slot Y_committee a list describing the models selected for
#' predictions at each time point (according to \strong{lambda}
#' and \strong{omega});
#'
#' @slot W a matrix with the weights of the base models at each prediction
#' time.
#'
#' @family ensemble predictions
#'
#' @export
setClass("dets_hat",
         slots = c(y_hat = "numeric",
                   Y_hat = "data.frame",
                   Y_committee = "OptionalList",
                   W = "data.frame")
)


#' Predictions by an DETS ensemble
#'
#' @param y_hat combined predictions of the ensemble
#' \code{\link{DETS}}. A numeric vector;
#'
#' @param Y_hat a matrix containing the predictions made by
#' individual models;
#'
#' @param Y_committee a list describing the models selected for
#' predictions at each time point (according to \strong{lambda}
#' and \strong{omega});
#'
#' @param W a matrix with the weights of the base models at each prediction
#' time.
#'
#' @family ensemble predictions
#'
#' @return Set of results from predicting with a \code{DETS}
#' ensemble
#'
#' @export
dets_hat <- function(y_hat, Y_hat, Y_committee, W) {

  if (missing(Y_committee)) Y_committee <- NULL

  methods::new(
    "dets_hat",
    y_hat = y_hat,
    Y_hat = Y_hat,
    Y_committee = Y_committee,
    W = W
  )
}

#' Forecasting using an ensemble predictive model
#'
#' Generic function for forecasting
#' future values of a time series from an \code{\link{ADE-class}} model or a
#' \code{\link{DETS-class}} model.
#'
#' @seealso \code{\link{predict}} for the predict method; \code{\link{update_weights}}
#' for updating the weights of a model after forecasting; \code{\link{update_base_models}}
#' for updating the base models of an ensemble.
#'
#' @param object predictive model object. A \code{\link{ADE-class}}
#' or a \code{\link{DETS-class}} ensemble object;
#'
#' @param h steps to forecast
#'
#' @family forecasting using ensembles
#'
#' @note the \code{forecast} generic in \strong{tsensembler} assumes that the
#' data is purely auto-regressive (no external variables), and that the target variable
#' is the first column of the data provided. For a different data setup, the
#' predict methods (\code{\link{predict}})
#' can be used (with successive calls with updates for multi-step forecasting).
#'
#' @examples
#' \dontrun{
#' specs <- model_specs(
#'  learner = c("bm_svr", "bm_glm", "bm_mars"),
#'  learner_pars = NULL
#' )
#'
#' data("water_consumption")
#' dataset <- embed_timeseries(water_consumption, 5)
#' train <- dataset[1:1000, ]
#'
#' model <- DETS(target ~., train, specs)
#' model2 <- ADE(target ~., train, specs, lambda = 30)
#'
#' next_vals_dets <- forecast(model, h = 2)
#' next_vals_ade <- forecast(model2, h = 2)
#' }
#'
#' @export
setGeneric("forecast",
           function(object, h) {
             standardGeneric("forecast")
           })
