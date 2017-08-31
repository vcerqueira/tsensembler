#' @title Predicting new observations using an ensemble
#'
#' @name predict
#' @aliases predict.ade predict.dets predict.base
#' @rdname predict-methods
#'
#' @description Initially, the predictions of the base models are collected.
#' Then, the predictions of the loss to be incurred by the base models \strong{E_hat}
#' (estimated by their associate meta models) are computed. The weights of
#' the base models are then estimated according to \strong{E_hat} and the committee of
#' top models. The committee is built according to the \emph{lambda} and
#' \emph{omega} parameters. Finally, the predictions are combined
#' according to the weights and the committee setup.
#'
#' @param object an object of class \code{\link{ADE-class}};
#' @param newdata new data to predict
#'
#'
#' @export
NULL

#' @rdname predict-methods
#'
#' @examples
#'
#' ###### Predicting with an ADE ensemble
#'
#' specs <- model_specs(
#'  learner = c("bm_glm", "bm_mars"),
#'  learner_pars = NULL
#' )
#'
#' data("water_consumption")
#' dataset <- embed_timeseries(water_consumption, 5)
#' train <- dataset[1:1000, ]
#' test <- dataset[1001:1500, ]
#'
#' model <- ADE(target ~., train, specs)
#'
#' preds <- predict(model, test)
#'
#'
#' @export
setMethod("predict",
          signature("ADE"),
          function(object, newdata) {
            seq. <- seq_len(nrow(newdata))

            Y_hat <- predict(object@base_ensemble, newdata)
            Y <- get_y(newdata, object@form)

            E_hat <- lapply(object@meta_model,
                            function(o)
                              predict(o, newdata)$predictions)

            E_hat <- as.data.frame(E_hat)

            W <- apply(E_hat,
                       1,
                       model_weighting,
                       trans = object@aggregation,
                       na.rm = TRUE)
            W <- t(W)

            if (!(object@select_best || object@all_models)) {
              Y_hat_recent <- predict(object@base_ensemble, object@recent_series)
              Y_rs <- get_y(object@recent_series, object@form)

              C <-
                build_committee(
                  rbind.data.frame(Y_hat_recent, Y_hat),
                  c(Y_rs, Y),
                  lambda = object@lambda,
                  omega = object@omega
                )

              C <- C[-seq_len(object@lambda)]
            } else
              C <- NULL

            if (object@select_best) {
              W <- select_best(W)
              C <- NULL
            }

            y_hat <-
              combine_predictions(
                Y_hat = Y_hat,
                W = W,
                committee = C
              )

            ade_hat(y_hat, Y_hat, C, E_hat)
          })


#' @rdname predict-methods
#'
#' @examples
#' \dontrun{
#'
#' ###### Predicting with a DETS ensemble
#'
#' specs <- model_specs(
#'  learner = c("bm_svr", "bm_glm", "bm_mars"),
#'  learner_pars = NULL
#' )
#'
#' data("water_consumption")
#' dataset <- embed_timeseries(water_consumption, 5)
#' train <- dataset[1:700, ]
#' test <- dataset[701:1000, ]
#'
#' model <- DETS(target ~., train, specs, lambda = 50, omega = .2)
#'
#' preds <- predict(model, test)
#' }
#'
#' @export
setMethod("predict",
          signature("DETS"),
          function(object, newdata) {
            seq. <- seq_len(nrow(newdata))

            Y_hat <- predict(object@base_ensemble, newdata)
            Y <- get_y(newdata, object@form)

            Y_hat_recent <-
              predict(object@base_ensemble, object@recent_series)
            Y_recent <- get_y(object@recent_series, object@form)

            scores <-
              model_recent_performance(rbind.data.frame(Y_hat_recent, Y_hat),
                                       c(Y_recent, Y),
                                       object@lambda,
                                       object@omega,
                                       object@base_ensemble@pre_weights)

            scores$top_models <- scores$top_models[-seq_len(object@lambda)]
            scores$model_scores <- scores$model_scores[-seq_len(object@lambda), ]

            y_hat <-
              combine_predictions(
                Y_hat = Y_hat,
                W = scores$model_scores,
                committee = scores$top_models
              )

            dets_hat(
              y_hat = y_hat,
              Y_hat = Y_hat,
              Y_committee = scores$top_models,
              W = scores$model_scores
            )
          })


#' @rdname predict-methods
#'
#' @examples
#'
#' \dontrun{
#' ###### Predicting with a base ensemble
#'
#' model <- ADE(target ~., train, specs)
#'
#' basepreds <- predict(model@base_ensemble, test)
#' }
#'
#'
#' @export
setMethod("predict",
          signature("base_ensemble"),
          function(object, newdata) {
            seq. <- seq_len(object@N)

            Y_hat <-
              lapply(seq.,
                     function(o) {
                       compute_predictions(object@base_models[o],
                                           object@form,
                                           newdata)
                     })

            Y_hat <- unlist(Y_hat, recursive = FALSE)
            nmodels <- names(Y_hat)

            Y_hat <- as.data.frame(Y_hat)
            colnames(Y_hat) <- nmodels

            Y_hat
          })
