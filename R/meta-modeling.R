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
           "ffnn" = {
             meta_ffnn(form, data)
           },
           "mars" = {
             meta_mars(form, data)
           },
           "ppr" = {
             meta_ppr(form, data)
           },
           "svr" = {
             meta_svr(form, data)
           },
           "pls" = {
             meta_pls(form, data)
           },
           "lasso" = {
             meta_lasso(form, data)
           },
           "gaussianprocess" = {
             meta_gp(form, data)
           },
           "xgb" = {
             meta_xgb(form, data)
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
           "ffnn" = {
             meta_ffnn_predict(model, newdata)
           },
           "svr" = {
             meta_svr_predict(model, newdata)
           },
           "mars" = {
             meta_mars_predict(model, newdata)
           },
           "pls" = {
             meta_pls_predict(model, newdata)
           },
           "ppr" = {
             meta_ppr_predict(model, newdata)
           },
           "lasso" = {
             meta_lasso_predict(model, newdata)
           },
           "gaussianprocess" = {
             meta_gp_predict(model, newdata)
           },
           "xgb" = {
             meta_xgb_predict(model, newdata)
           },
           stop("unknown meta-model.")
    )
  }


####################################################################
####################################################################

#' Training a Gaussian process arbiter
#'
#' @param form form
#' @param data data
#'
#' @keywords internal
#'
#' @export
meta_svr <-
  function(form, data) {
    bm_svr(form, data, NULL)[[1]]
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
meta_svr_predict <-
  function(model, newdata) {
    predict(model, newdata)[,1]
  }


#' Training a pls process arbiter
#'
#' @param form form
#' @param data data
#'
#'
#' @keywords internal
#'
#' @export
meta_pls <-
  function(form, data) {
    model <- pls::mvr(formula = form,
                      data = data,
                      method = "kernelpls")

    model$best_comp_train <- best_mvr(model, form, data)

    model
  }

#' Arbiter predictions via pls model
#'
#' @param model arbiter -- a Gaussian process model
#' @param newdata new data to predict loss
#'
#' @keywords internal
#'
#' @export
meta_pls_predict <-
  function(model, newdata) {
    predict_pls_pcr(model, newdata)
  }


#' Training a meta_mars process arbiter
#'
#' @param form form
#' @param data data
#'
#' @keywords internal
#'
#' @export
meta_ppr <-
  function(form, data) {
    stats::ppr(form,
               data,
               nterms = 5,
               sm.method = "supsmu")
  }

#' Arbiter predictions via ppr model
#'
#' @param model arbiter -- a Gaussian process model
#' @param newdata new data to predict loss
#'
#' @keywords internal
#'
#' @export
meta_ppr_predict <-
  function(model, newdata) {
    unname(predict(model, newdata))
  }



#' Training a meta_mars process arbiter
#'
#' @param form form
#' @param data data
#'
#' @import kernlab
#'
#' @keywords internal
#'
#' @export
meta_mars <-
  function(form, data) {
    earth::earth(form,
                 data,
                 nk = 10,
                 degree = 3,
                 thresh = 0.001)
  }

#' Arbiter predictions via mars model
#'
#' @param model arbiter -- a Gaussian process model
#' @param newdata new data to predict loss
#'
#' @import earth
#'
#' @keywords internal
#'
#' @export
meta_mars_predict <-
  function(model, newdata) {
    predict(model, newdata)[,1]
  }


#' Training a Gaussian prosadacess arbiter
#'
#' @param form form
#' @param data data
#'
#' @keywords internal
#'
#' @export
meta_ffnn <-
  function(form, data) {
    tr_X <- stats::model.matrix(form, data)
    tr_Y <- get_y(data, form)

    monmlp::monmlp.fit(tr_X,
                       t(t(tr_Y)),
                       hidden1=10,
                       hidden2=0,
                       n.ensemble=1,
                       bag=F,
                       silent=T)
  }

#' Arbiter predictions via linear ssmodel
#'
#' @param model arbiter -- a Gaussian process model
#' @param newdata new data to predict loss
#'
#'
#' @keywords internal
#'
#' @export
meta_ffnn_predict <-
  function(model, newdata) {
    form <- score ~.

    ts_X <- as.matrix(stats::model.matrix(form, newdata))
    preds <- monmlp::monmlp.predict(ts_X, model)[,1]

    unname(preds)
  }



##



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
            kernel = "vanilladot",
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


###

#' Training a xgb arbiter
#'
#' @param form formula
#' @param data data
#'
#' @keywords internal
#'
#' @export
meta_xgb <-
  function(form, data) {
    model <- bm_xgb(form, data)
    model <- model$xgb_optimd

    model
  }

#' Arbiter predictions via xgb
#'
#' @param meta_model arbiter -- a ranger object
#' @param newdata new data to predict
#'
#' @export
meta_xgb_predict <-
  function(meta_model, newdata) {
    xgb_predict_(meta_model, newdata)
  }

#' asdasd
#'
#' @param model mode
#' @param newdata s
#'
#' @import xgboost
#'
#' @export
xgb_predict_ <- function(model, newdata) {
  id <- which(colnames(newdata) == "target")
  colnames(newdata)[id] <- "score"

  Y <- get_y(newdata, model$form)

  valid_cols <- which(colnames(newdata) %in% model$cn)

  newdata  <- subset(newdata, select = valid_cols)
  newdata  <- stats::model.matrix(model$form, newdata)

  dtest <- xgb.DMatrix(newdata, label = Y)

  predict(model, dtest)
}

###

