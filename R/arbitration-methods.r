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
#' @param lfun meta loss function
#'
#' @keywords internal
#'
#' @export
"train_ade" <-
  function(form, train, specs, lambda, lfun) {
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
        lfun = lfun
      )

    cat("Building base ensemble...\n")
    M <- build_base_ensemble(form, train, specs)

    OOB <- rbind_l(lapply(OOB_data, function(x) x$oob))
    train_loss <- rbind_l(lapply(OOB_data, function(x) x$mloss))
    recent_lambda_k <- recent_lambda_observations(OOB, lambda)

    OOB_no_tgt <- subset(OOB, select = -which(colnames(OOB) %in% tgt))

    train_metadata <-
      lapply(train_loss,
             function(l_hat) cbind.data.frame(OOB_no_tgt,
                                              score = l_hat))

    cat("Training meta models (an arbiter for each expert\n")
    meta_models <-
      lapply(train_metadata,
             function(meta_set) {
              if (any(is.na(meta_set))) {
                meta_set <- soft.completion(meta_set)
              }
               loss_meta_learn(score ~ ., meta_set, "randomforest")
             })

    list(
      base_ensemble = M,
      meta_model = meta_models,
      recent_series = recent_lambda_k,
      OOB = OOB_data
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

            unique_M <- unique(vcapply(names_M, function(j) split_by_(j)[1]))
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


#' @rdname forecast
setMethod("forecast",
          signature("ADE"),
          function(object, h) {

            forecasts <- numeric(h)
            for (o in seq_len(h)) {
              p_forecast <- forecast_ade(object)
              forecasts[o] <- p_forecast

              object@recent_series <-
                rbind.data.frame(object@recent_series,
                                 struc_embed(object, p_forecast))
            }
            names(forecasts) <- paste("t", seq_along(forecasts), sep = "+")

            forecasts
          })

forecast_ade <-
  function(object) {
    next_embed <- struc_embed(object, -1)

    Y_hat <- predict(object@base_ensemble, next_embed)

    E_hat <- lapply(object@meta_model,
                    function(o)
                      predict(o, next_embed)$predictions)

    E_hat <- as.data.frame(E_hat)

    W <- t(apply(E_hat, 1, model_weighting, na.rm = TRUE))

    Y_hat_recent <- predict(object@base_ensemble, object@recent_series)

    Y_rs <- get_y(object@recent_series, object@form)

    C <- build_committee(Y_hat_recent, Y_rs, object@lambda, object@omega)
    last_C <- C[[length(C)]]

    Y_hat_j <- Y_hat[1, last_C]
    W_j <- proportion(W[1, last_C])

    point_forecast <- sum(Y_hat_j * W_j)

    point_forecast
  }


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
           function(object, newdata) {
             standardGeneric("update_ade_meta")
           })

#' @rdname update_ade_meta
setMethod("update_ade_meta",
          signature("ADE"),
          function(object, newdata) {
            OOB_data <-
              blocked_prequential(
                x = newdata,
                nfolds = 10,
                intraining_estimations,
                .rbind = FALSE,
                form = object@form,
                specs = object@specs,
                lfun = ae
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
           function(object, newdata) {
             standardGeneric("update_ade")
           })

#' @rdname update_ade
setMethod("update_ade",
          signature("ADE"),
          function(object, newdata) {

            object <- update_base_models(object, newdata)
            object <- update_ade_meta(object, newdata)

            object
          })
