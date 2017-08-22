#' Updating the weights of base models
#'
#' Update the weights of base models of a \code{\link{ADE-class}}
#' or \code{\link{DETS-class}} ensemble.
#' This is accomplished by using computing the loss of the base models
#' in new recent observations.
#'
#' @param object a \code{\link{ADE-class}} or \code{\link{DETS-class}} model object;
#'
#' @param newdata new data used to update the most
#' recent observations of the time series. At prediction time
#' these observations are used to compute the weights of the base models
#'
#' @note Updating the weights of an ensemble is only necessary between
#' different calls of the functions \code{predict} or \code{forecast}.
#' Otherwise, if consecutive know observations are predicted
#' (e.g. a validation/test set) the updating is automatically done internally.
#'
#' @family updating models
#'
#' @seealso \code{\link{update_weights}} for the weight updating method
#' for an \code{\link{ADE}} model, and \code{\link{update_weights}} for the same method
#' for a \code{\link{DETS}} model
#'
#' @examples
#' \dontrun{
#' data("water_consumption")
#' dataset <- embed_timeseries(water_consumption, 5)
#'
#' train <- dataset[1:1000,]
#' test <- dataset[1001:1010, ]
#'
#' specs <- model_specs(c("bm_ppr","bm_glm","bm_mars"), NULL)
#' ## same with model <- DETS(target ~., train, specs)
#' model <- ADE(target ~., train, specs)
#'
#' # if consecutive know observations are predicted (e.g. a validation/test set)
#' # the updating is automatically done internally.
#' predictions1 <- predict(model, test)@y_hat
#'
#' # otherwise, the models need to be updated
#' predictions <- numeric(nrow(test))
#' # predict new data and update the weights of the model
#' for (i in seq_along(predictions)) {
#'   predictions[i] <- predict(model, test[i, ])@y_hat
#'
#'   model <- update_weights(model, test[i, ])
#' }
#'
#' all.equal(predictions1, predictions)
#' }
#'
#' @export
setGeneric("update_weights",
           function(object, newdata) {
             standardGeneric("update_weights")
           })

#' @rdname update_weights
setMethod("update_weights",
          signature("ADE"),
          function(object, newdata) {
            rseries <- object@recent_series

            if (!all(colnames(object@recent_series) == colnames(newdata)))
              stop(
                "The dimension of the new data must have the
                same structure as the \"object@recent_series\" (colnames)",
                call. = FALSE
              )

            rseries <- rbind.data.frame(rseries, newdata)

            rseries <- recent_lambda_observations(rseries, object@lambda)

            object@recent_series <- rseries

            object
          })

#' @rdname update_weights
setMethod("update_weights",
          signature("DETS"),
          function(object, newdata) {
            rseries <- object@recent_series

            if (!all(colnames(object@recent_series) == colnames(newdata)))
              stop(
                "The dimension of the new data must have the
                same structure as the \"object@recent_series\" (colnames)",
                call. = FALSE
              )

            rseries <- rbind.data.frame(rseries, newdata)

            rseries <- recent_lambda_observations(rseries, object@lambda)

            object@recent_series <- rseries

            object
          })
