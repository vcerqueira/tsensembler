#' xgb predict fun
#' @import xgboost
#'
#' @param model model from bm_xgb
#' @param newdata newdata
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

    dtest <- xgb.DMatrix(newdata, label = Y)

    predict(model, dtest)
  }


#' xgb base model
#'
#' @param form formula
#' @param data training data
#' @param lpars learning parameters
#'
#' @keywords internal
#'
#' @import xgboost
#' @export
bm_xgb <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_xgb))
      lpars$bm_xgb <- list()

    if (is.null(lpars$bm_xgb$max_depth))
      lpars$bm_xgb$max_depth <- 4

    if (is.null(lpars$bm_xgb$eta))
      lpars$bm_xgb$eta <- .1

    if (is.null(lpars$bm_xgb$nrounds))
      lpars$bm_xgb$nrounds <- 100

    nmodels <-
      length(lpars$bm_xgb$max_depth) *
      length(lpars$bm_xgb$eta) *
      length(lpars$bm_xgb$nrounds)

    cn <- colnames(data)
    Y <- get_y(data, form)
    Y <- as.integer(as.character(Y))
    data <- model.matrix.na(form, data)
    dtrain <- xgb.DMatrix(data, label = Y)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (max_depth in lpars$bm_xgb$max_depth) {
      for (eta in lpars$bm_xgb$eta) {
        for (nrounds in lpars$bm_xgb$nrounds) {
          j <- j + 1L

          fparams <-
            list(
              max_depth = max_depth,
              eta = eta,
              silent = 1,
              objective = "reg:linear",
              eval_metric = "rmse"
            )

          ensemble[[j]] <-
            xgb.train(fparams,
                      dtrain,
                      nrounds = nrounds)

          #ensemble[[j]] <-
          #  xgb.Booster.complete(ensemble[[j]])

          ensemble[[j]]$cn <- cn
          ensemble[[j]]$form <- form

          #ensemble[[j]] <- list(ensemble[[j]])

          paramsnames <- paste0("md_", max_depth,
                                "e_", eta,
                                "nr_", nrounds)

          mnames[[j]] <- paste0("xgb_", paramsnames)
        }
      }
    }

    names(ensemble) <- mnames
    #cat("here\n")

    ensemble
  }
