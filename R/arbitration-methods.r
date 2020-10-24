#' Training procedure of for ADE
#'
#' Base level models are trained according to \strong{specs},
#' and meta level models are trained using a blocked prequential
#' procedure in out-of-bag samples from the training data.
#'
#' @param form formula;
#'
#' @param train training data as a data frame;
#'
#' @param specs a \code{\link{model_specs-class}} object class. It contains
#' the parameter setting specifications for training the ensemble;
#'
#' @param lambda window size. Number of observations to compute
#' the recent performance of the base models, according to the
#' committee ratio \strong{omega}. Essentially, the top \emph{omega}
#' models are selected and weighted at each prediction instance, according
#' to their performance in the last \emph{lambda} observations.
#' Defaults to 50 according to empirical experiments;
#'
#' @param lfun meta loss function - defaults to \strong{ae} (absolute error)
#'
#' @param meta_model_type algorithm used to train meta models. Defaults to a
#' random forest (using ranger package)
#'
#' @param num_cores A numeric value to specify the number of cores used to
#' train base and meta models. num_cores = 1
#' leads to sequential training of models. num_cores > 1
#' splits the training of the base models across num_cores cores.
#'
#' @keywords internal
#'
#' @import parallel
#' @import foreach
#'
#' @export
"train_ade" <-
  function(form, train, specs, lambda, lfun, meta_model_type, num_cores) {
    if (length(num_cores) > 1 || !is.numeric(num_cores)) {
      stop("Please specify a numeric value for num_cores. num_cores = 1
           leads to sequential training of models. num_cores > 1
           splits the training of the base models across num_cores cores.")
    }
    if (is.null(num_cores)) num_cores <- 1

    tgt <- get_target(form)

    cat("Setting up meta data\n")
    OOB_data <-
      blocked_prequential(
        x = train,
        nfolds = 10,
        intraining_estimations,
        .rbind = FALSE,
        form = form,
        specs = specs,
        lfun = lfun,
        num_cores = num_cores
      )

    cat("Building base ensemble...\n")
    M <- build_base_ensemble(form, train, specs, num_cores)

    OOB <- rbind_l(lapply(OOB_data, function(x) x$oob))
    train_loss <- rbind_l(lapply(OOB_data, function(x) x$mloss))
    recent_lambda_k <- recent_lambda_observations(OOB, lambda)

    OOB_no_tgt <- subset(OOB, select = -which(colnames(OOB) %in% tgt))

    train_metadata <-
      lapply(train_loss,
             function(l_hat) cbind.data.frame(OOB_no_tgt,
                                              score = l_hat))

    cat("Training meta models (an arbiter for each expert\n")
    if (num_cores > 1) {
      cl <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl)
      `%partrain%` <- `%dopar%`
    } else {
      `%partrain%` <- `%do%`
    }

    meta_models <-
      foreach::foreach(meta_set = train_metadata,
                       .packages = "tsensembler") %partrain% {

                         if (any(is.na(meta_set))) {
                           meta_set <- soft.completion(meta_set)
                         }

                         loss_meta_learn(score ~ ., meta_set, meta_model_type)
                       }

    list(
      base_ensemble = M,
      meta_model = meta_models,
      recent_series = recent_lambda_k,
      OOB = OOB_data
    )
  }

#' ADE training poor version
#' Train meta-models in the training data,
#' as opposed to using a validation dataset
#'
#' Saves times by not computing oob predictions.
#' Testing comp costs are the same.
#'
#' @param form formula
#' @param train training data
#' @param specs a \code{\link{model_specs-class}} object class. It contains
#' the parameter setting specifications for training the ensemble;
#'
#' @param lambda window size. Number of observations to compute
#' the recent performance of the base models, according to the
#' committee ratio \strong{omega}. Essentially, the top \emph{omega}
#' models are selected and weighted at each prediction instance, according
#' to their performance in the last \emph{lambda} observations.
#' Defaults to 50 according to empirical experiments;
#'
#' @param lfun meta loss function - defaults to \strong{ae} (absolute error)
#'
#' @param meta_model_type algorithm used to train meta models. Defaults to a
#' random forest (using ranger package)
#'
#' @param num_cores A numeric value to specify the number of cores used to
#' train base and meta models. num_cores = 1
#' leads to sequential training of models. num_cores > 1
#' splits the training of the base models across num_cores cores.
#'
#' @keywords internal
#'
#' @import parallel
#' @import foreach
#'
#' @export
train_ade_quick <-
  function(form, train, specs, lambda, lfun, meta_model_type, num_cores) {
    if (length(num_cores) > 1 || !is.numeric(num_cores)) {
      stop("Please specify a numeric value for num_cores. num_cores = 1
           leads to sequential training of models. num_cores > 1
           splits the training of the base models across num_cores cores.")
    }
    if (is.null(num_cores)) num_cores <- 1

    tgt <- get_target(form)

    cat("Building base ensemble...\n")
    M <- build_base_ensemble(form, train, specs, num_cores)

    cat("Setting up poor meta data \n")
    y_tr <- get_y(train, form)
    Y_hat_tr <- predict(M, train)

    train_loss <-
      base_models_loss(
        Y_hat = Y_hat_tr,
        Y = y_tr,
        lfun = lfun)

    OOB <- train
    OOB_no_tgt <- subset(OOB, select = -which(colnames(OOB) %in% tgt))

    recent_lambda_k <- recent_lambda_observations(OOB, lambda)

    train_metadata <-
      lapply(train_loss,
             function(l_hat) cbind.data.frame(OOB_no_tgt,
                                              score = l_hat))

    cat("Training meta models (an arbiter for each expert\n")
    if (num_cores > 1) {
      cl <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl)
      `%partrain%` <- `%dopar%`
    } else {
      `%partrain%` <- `%do%`
    }

    meta_models <-
      foreach::foreach(meta_set = train_metadata,
                       .packages = "tsensembler") %partrain% {

                         if (any(is.na(meta_set))) {
                           meta_set <- soft.completion(meta_set)
                         }

                         loss_meta_learn(score ~ ., meta_set, meta_model_type)
                       }

    list(
      base_ensemble = M,
      meta_model = meta_models,
      recent_series = recent_lambda_k,
      OOB = OOB
    )
    }


setMethod("show",
          signature("ADE"),
          function(object) {
            cat("## Arbitrated Dynamic Ensemble ##\n\n")

            cat("Base ensemble specs:\n")
            print(object@base_ensemble)
            cat("\n\n")

            cat("Lambda is set to:", object@lambda, "\n")
            if (object@select_best)
              cat("Selecting the predicted best base model at each point.\n")
            else if (object@all_models)
              cat("Combining all models, using a ", object@aggregation, "function.")
            else
              cat("Omega is set to:", object@omega, "\n\n")


          }
)

setMethod("show",
          signature("ade_hat"),
          function(object) {
            cat("Predictions from an ADE model\n\n")
            names_M <- names(object@Y_hat)
            len <- length(names_M)

            unique_M <- vapply(names_M, function(j) {
                split_by_(j)[1]
              }, character(1), USE.NAMES=FALSE)

            unique_M <- unique(unique_M)
            unique_M <- paste(unique_M, collapse = ", ")

            cat("The predictions were calculated with", len, "base models:\n")
            cat(unique_M, "\n\n")
            cat("Check:\n")
            cat("slot @y_hat for the combined predictions\n")
            cat("slot @Y_hat for the predictions of each base model\n")
            cat("slot @Y_committee to see the selected models at each prediction\n")
            cat("slot @E_hat for the predictions of each meta model\n")
          }
)

#' Updating the metalearning layer of an ADE model
#'
#' The \strong{update_ade_meta} function uses new information to
#' update the meta models of an \code{\link{ADE-class}} ensemble. As input
#' it receives a \code{\link{ADE-class}} model object class and a new dataset
#' for updating the weights of the base models in the ensemble.
#' This new data should have the same structure as the one used to build the
#' ensemble. Updating the base models of the ensemble is done using the \code{\link{update_base_models}}
#' function.
#'
#' @param object a \code{\link{ADE-class}} object.
#'
#' @param newdata data used to update the meta models. This should be
#' the data used to initially train the meta-models (training set), together
#' with new observations (for example, validation set). Each meta model
#' is retrained using \code{newdata}.
#'
#' @param num_cores A numeric value to specify the number of cores used to
#' train base and meta models. num_cores = 1
#' leads to sequential training of models. num_cores > 1
#' splits the training of the base models across num_cores cores.
#'
#' @family updating models
#'
#' @seealso \code{\link{ADE-class}} for building an ADE model;
#' \code{\link{update_weights}} for updating the weights of the ensemble (without
#' retraining the models); and \code{\link{update_base_models}} for updating the
#' base models of an ensemble.
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
#' validation <- dataset[1001:1200, ]
#' test <- dataset[1201:1500, ]
#'
#' model <- ADE(target ~., train, specs)
#'
#' preds_val <- predict(model, validation)
#' model <- update_ade_meta(model, rbind.data.frame(train, validation))
#'
#' preds_test <- predict(model, test)
#' }
#'
#' @export
setGeneric("update_ade_meta",
           function(object, newdata, num_cores=1) {
             standardGeneric("update_ade_meta")
           })

#' @rdname update_ade_meta
setMethod("update_ade_meta",
          signature("ADE"),
          function(object, newdata, num_cores=1) {
            OOB_data <-
              blocked_prequential(
                x = newdata,
                nfolds = 10,
                intraining_estimations,
                .rbind = FALSE,
                form = object@form,
                specs = object@specs,
                lfun = ae,
                num_cores=num_cores
              )

            tgt <- get_target(object@form)

            OOB <- rbind_l(lapply(OOB_data, function(x) x$oob))
            recent_lambda_k <-
              recent_lambda_observations(OOB, object@lambda)

            OOB <- subset(OOB, select = -which(colnames(OOB) %in% tgt))

            train_loss <- rbind_l(lapply(OOB_data, function(x) x$mloss))

            train_metadata <-
              lapply(train_loss,
                     function(l_hat) cbind.data.frame(OOB,
                                                      score = l_hat))

            meta_models <-
              lapply(train_metadata,
                     function(meta_set) {
                      if (any(is.na(meta_set))) {
                        meta_set <- soft.completion(meta_set)
                      }

                      ranger(score ~ .,
                             meta_set,
                             mtry = NCOL(meta_set) / 3,
                             num.trees = 500,
                             write.forest = TRUE)
                     })

            object@meta_model <- meta_models
            object@recent_series <- recent_lambda_k

            object
          })


#' Updating an ADE model
#'
#' \strong{update_ade} is a generic function that combines
#' \code{\link{update_base_models}}, \code{\link{update_ade_meta}},
#' and \code{\link{update_weights}}.
#'
#' @param object a \code{\link{ADE-class}} object.
#'
#' @param newdata data used to update the ADE model. This should be
#' the data used to initially train the models (training set), together
#' with new observations (for example, validation set). Each model
#' is retrained using \code{newdata}.
#'
#' @param num_cores A numeric value to specify the number of cores used to
#' train base and meta models. num_cores = 1
#' leads to sequential training of models. num_cores > 1
#' splits the training of the base models across num_cores cores.
#'
#' @family updating models
#'
#' @seealso \code{\link{ADE-class}} for building an ADE model;
#' \code{\link{update_weights}} for updating the weights of the ensemble (without
#' retraining the models); \code{\link{update_base_models}} for updating the
#' base models of an ensemble; and \code{\link{update_ade_meta}} for
#' updating the meta-models of an ADE model.
#'
#' @examples
#' specs <- model_specs(
#'  learner = c("bm_svr", "bm_glm", "bm_mars"),
#'  learner_pars = NULL
#' )
#'
#' data("water_consumption")
#' dataset <- embed_timeseries(water_consumption, 5)
#' # toy size for checks
#' train <- dataset[1:300, ]
#' validation <- dataset[301:400, ]
#' test <- dataset[401:500, ]
#'
#' model <- ADE(target ~., train, specs)
#'
#' preds_val <- predict(model, validation)
#' model <- update_ade(model, rbind.data.frame(train, validation))
#'
#' preds_test <- predict(model, test)
#'
#'
#' @export
setGeneric("update_ade",
           function(object, newdata, num_cores=1) {
             standardGeneric("update_ade")
           })

#' @rdname update_ade
setMethod("update_ade",
          signature("ADE"),
          function(object, newdata,num_cores=1) {

            object <- update_base_models(object, newdata,num_cores)
            object <- update_ade_meta(object, newdata,num_cores)

            object
          })
