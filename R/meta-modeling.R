#' Arbiter predictions via ranger
#'
#' @param meta_model arbiter -- a ranger object
#' @param newdata new data to predict
#'
#' @import ranger
#'
#' @keywords internal
#'
#' @export
meta_rf_predict <-
  function(meta_model, newdata) {
    predict(meta_model, newdata)$predictions
  }

#' Arbiter predictions via linear model
#'
#' @param meta_model arbiter -- a glmnet object
#' @param newdata new data to predict
#'
#' @keywords internal
#'
#' @export
meta_lasso_predict <-
  function(meta_model, newdata) {
    tgt <- get_target(meta_model$form)

    X_names <- meta_model$beta@Dimnames[[1]]
    newdata <- newdata[colnames(newdata) %in% c(tgt, X_names)]
    newdata <- cbind.data.frame(score = -1., newdata)

    X <- stats::model.matrix(meta_model$form, newdata)

    predict.glmnet(meta_model, X, type = "response")[,1]
  }

#' Training a random forest arbiter
#'
#' @param form formula
#' @param data data
#'
#' @import ranger
#'
#' @keywords internal
#'
#' @export
meta_rf <-
  function(form, data) {
    ranger(form,
           data,
           mtry = NCOL(data) / 3,
           num.trees = 500,
           write.forest = TRUE,
           keep.inbag = TRUE)
  }

#' Training a LASSO arbiter
#'
#' @param form form
#' @param data data
#'
#' @import glmnet
#'
#' @keywords internal
#'
#' @export
meta_lasso <-
  function(form, data) {
    alpha <- 1

    X <- stats::model.matrix(form, data)
    Y <- get_y(data, form)

    m.all <- glmnet(X, Y, alpha = alpha)

    model <- glmnet(X, Y, alpha = 1, lambda = min(m.all$lambda))

    model$form <- form

    model
  }


#' Training an arbiter
#'
#' @param form form
#' @param data data
#' @param meta_model learning algorithm -- either a "randomforest",
#' a "lasso", or a "gaussianprocess".
#'
#' @keywords internal
#'
#'
#' @export
loss_meta_learn <-
  function(form, data, meta_model) {
    switch(meta_model,
           "randomforest" = {
             meta_rf(form, data)
           },
           "cubist" = {
             meta_cubist(form, data)
           },
           "lasso" = {
             meta_lasso(form, data)
           },
           "gaussianprocess" = {
             meta_gp(form, data)
            },
           stop("unknown meta-model.")
           )
  }

#' Predicting loss using arbiter
#'
#' @param model arbiter model
#' @param newdata new data to predict loss
#' @param meta_model learning algorithm -- either a "randomforest",
#' a "lasso", or a "gaussianprocess".
#'
#'
#' @keywords internal
#'
#' @export
meta_predict <-
  function(model, newdata, meta_model) {
    switch(meta_model,
           "randomforest" = {
             meta_rf_predict(model, newdata)
           },
           "cubist" = {
             meta_cubist_predict(model, newdata)
           },
           "lasso" = {
             meta_lasso_predict(model, newdata)
           },
           "gaussianprocess" = {
             meta_gp_predict(model, newdata)
           },
           stop("unknown meta-model.")
    )
  }

#' Training a Gaussian process arbiter
#'
#' @param form form
#' @param data data
#'
#' @import kernlab
#'
#' @keywords internal
#'
#' @export
meta_gp <-
  function(form, data) {
    kernlab::gausspr(form,
            data,
            type = "regression",
            kernel = "rbfdot",
            tol = .01)
  }

#' Arbiter predictions via linear model
#'
#' @param model arbiter -- a Gaussian process model
#' @param newdata new data to predict loss
#'
#'
#' @import kernlab
#'
#' @keywords internal
#'
#' @export
meta_gp_predict <-
  function(model, newdata) {
    predict(model, newdata)[,1]
  }


#' Training a RBR arbiter
#'
#' @param form formula
#' @param data data
#'
#' @import Cubist
#'
#' @keywords internal
#'
#' @export
meta_cubist <-
  function(form, data) {
    tr_X <- stats::model.matrix(form, data)
    tr_Y <- get_y(data, form)

    model <- suppressWarnings(Cubist::cubist(tr_X,tr_Y, committees = 50))
    model$form <- form

    model
  }

#' Arbiter predictions via Cubist
#'
#' @param meta_model arbiter -- a ranger object
#' @param newdata new data to predict
#'
#' @import Cubist
#'
#' @keywords internal
#'
#' @export
meta_cubist_predict <-
  function(meta_model, newdata) {

    vars <- meta_model$vars$all
    ids <- which(!colnames(newdata) %in% vars)

    newdata <- subset(newdata, select = -ids)
    newdata$score <- -1

    ts_X <- stats::model.matrix(meta_model$form, newdata)

    predict(meta_model, ts_X)
  }
