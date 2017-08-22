#' Update the base models of an ensemble
#'
#' This is a generic function for updating the base models
#' comprising an ensemble.
#'
#' \strong{update_base_models} function receives a
#' model object and a new dataset for retraining
#' the base models. This new data should
#' have the same structure as the one used to build the
#' ensemble.
#'
#' @param object an ensemble object, of class \code{\link{DETS-class}} or
#' \code{\link{ADE-class}};
#'
#' @param newdata new data used to update the models. Each base model
#' is retrained, so \code{newdata} should be the past data used for initially training
#' the models along with any further available observations.
#'
#' @seealso \code{\link{ADE-class}} for the ADE model information, and
#' \code{\link{DETS-class}} for the DETS model information;
#' \code{\link{update_ade_meta}} for updating the meta models of an ADE ensemble.
#' See \code{\link{update_weights}} for the method used to update
#' the weights of the ensemble. Updating the weights only changes the information
#' about the recent observations for computing the weights of the base models,
#' while updating the model uses that information to retrain the models.
#'
#' @examples
#' \dontrun{
#' data("water_consumption")
#' dataset <- embed_timeseries(water_consumption, 5)
#' train <- dataset[1:1000,]
#' test <- dataset[1001:1010, ]
#'
#' specs <- model_specs(c("bm_ppr","bm_glm","bm_mars"), NULL)
#'
#' model <- ADE(target ~., train, specs)
#'
#' predictions <- numeric(nrow(test))
#' for (i in seq_along(predictions)) {
#'   predictions[i] <- predict(model, test[i, ])@y_hat
#'   model <-
#'     update_base_models(model,
#'                        rbind.data.frame(train, test[seq_len(i), ]))
#' }
#'
#' ####
#'
#' specs2 <- model_specs(c("bm_ppr","bm_randomforest","bm_svr"), NULL)
#'
#' modeldets <- DETS(target ~., train, specs2)
#'
#' predictions <- numeric(nrow(test))
#' # predict new data and update models every three points
#' # in the remaining points, the only the weights are updated
#' for (i in seq_along(predictions)) {
#'   predictions[i] <- predict(modeldets, test[i, ])@y_hat
#'
#'   if (i %% 3 == 0)
#'     modeldets <-
#'       update_base_models(modeldets,
#'                          rbind.data.frame(train, test[seq_len(i), ]))
#'   else
#'     modeldets <- update_weights(modeldets, test[seq_len(i), ])
#' }
#' }
#'
#' @export
setGeneric("update_base_models",
           function(object, newdata) {
             standardGeneric("update_base_models")
           })

#' @rdname update_base_models
setMethod("update_base_models",
          signature("ADE"),
          function(object, newdata) {

            object@base_ensemble <-
              build_base_ensemble(object@form, newdata, object@specs)

            object <- update_weights(object, newdata)

            object
          })

#' @rdname update_base_models
setMethod("update_base_models",
          signature("DETS"),
          function(object, newdata) {
            recent_lambda_k <- recent_lambda_observations(newdata, object@lambda)

            object@base_ensemble <-
              build_base_ensemble(object@form, newdata, object@specs)

            object@recent_series <- recent_lambda_k

            object
          })
