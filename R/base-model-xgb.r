#' Base model for XGBoost
#'
#' @param form formula
#' @param data Training data
#' @param lpars list of parameters--deprecated
#'
#' @keywords internal
#'
#' @export
bm_xgb <-
  function(form, data, lpars) {
    GS <-
      expand.grid(
        eta = c(.1, .25, .5, .75),
        max_depth = c(1, 2, 3, 5, 10),
        booster = c("gblinear", "gbtree")
      )

    cn <- colnames(data)
    y <- get_y(data, form)
    X <- stats::model.matrix(form, data)

    Xy <- xgboost::xgb.DMatrix(X, label = y)

    params <- xgb_optimizer(X, y, GS)

    fparams <-
      list(
        max_depth = params$max_depth,
        eta = params$eta,
        silent = 1,
        nthread = 10,
        objective = "reg:squarederror",
        eval_metric = "rmse"
      )

    model <-
      xgboost::xgb.train(
        params = fparams,
        booster = as.character(params$booster),
        data = Xy,
        nrounds = params$nrounds,
        verbose = 0
      )

    model$cn <- cn
    model$form <- form

    model <- list(model)

    names(model) <- "xgb_optimd"

    model
  }


#' XGBoost predict function
#'
#' @param model Model from bm_xgb
#' @param newdata Test data
#'
#' @keywords internal
#'
#' @export
xgb_predict <-
  function(model, newdata) {
    Y <- get_y(newdata, model$form)
    valid_cols <- which(colnames(newdata) %in% model$cn)

    newdata  <- subset(newdata, select = valid_cols)
    newdata  <- stats::model.matrix(model$form, newdata)

    dtest <- xgboost::xgb.DMatrix(newdata, label = Y)

    predict(model, dtest)
  }



#' XGB optimizer
#'
#' @param X Covariates
#' @param y Target values
#' @param gsearch Grid search
#'
#' @export
xgb_optimizer <-
  function(X, y, gsearch) {
    l <- nrow(X)
    train_ids <- 1:ceiling(l * .7)

    X_tr <- X[train_ids, ]
    X_ts <- X[-train_ids, ]
    y_tr <- y[train_ids]
    y_ts <- y[-train_ids]

    tr <- xgboost::xgb.DMatrix(as.matrix(X_tr), label = y_tr)
    tst <- xgboost::xgb.DMatrix(as.matrix(X_ts), label = y_ts)

    wlist <- list(eval = tst, train = tr)

    L <- numeric(nrow(gsearch))
    NROUNDS <- numeric(nrow(gsearch))
    for (i in 1:nrow(gsearch)) {
      #i<-1
      cat(i, "\n")
      eta_i <- gsearch[i, "eta"]
      md_i <- gsearch[i, "max_depth"]
      algo_i <- as.character(gsearch[i, "booster"])

      params_i  <-
        list(
          max_depth = md_i,
          eta = eta_i,
          silent = 1,
          nthread = 1,
          objective = "reg:linear",
          eval_metric = "rmse"
        )

      model <-
        xgboost::xgb.train(
          params = params_i,
          booster = algo_i,
          data = tr,
          nrounds = 50,
          watchlist = wlist,
          verbose = 0
        )

      eval.log <- model$evaluation_log

      val.loss <- eval.log[[2]]
      L[i] <- min(val.loss)
      NROUNDS[i] <- which.min(val.loss)
    }

    best_var <- which.min(L)
    nrounds_best <- NROUNDS[best_var]

    eta <- gsearch[best_var, "eta"]
    max_depth <- gsearch[best_var, "max_depth"]
    booster <- as.character(gsearch[best_var, "booster"])

    best_parameters <-
      list(
        eta = eta,
        max_depth = max_depth,
        nrounds = nrounds_best,
        booster = booster
      )

    best_parameters
  }
