#' Fitting an ARIMA model
#' Using function \strong{auto.arima} from forecast package
#'
#' @param y time series struture after applying
#' function \strong{embed_timeseries}
#'
#' @param form formula
#'
#' @keywords internal
#'
#' @export
bm_arima <-
  function(y, form) {
    tm_ids <- c(1, grep("Tm", colnames(y)))
    train_xreg <- subset(y, select = -tm_ids)
    if (ncol(train_xreg) < 1) {
      train_xreg <- NULL
    } else {
      train_xreg <- as.matrix(train_xreg)
    }

    y <- get_y(y, form)
    model <- tryCatch(forecast::auto.arima(y = y, xreg = train_xreg),
                      error = function(e) NULL)
    if (is.null(model)) {
      model <- forecast::Arima(y, order = c(0,1,0))
    }
    model$form <- form

    model
  }

#' Predict function for a model created from bm_arima
#'
#' @param model model created using \strong{bm_arima} function
#'
#' @param test test data created using \strong{embed_timeseries} function
#'
#' @keywords internal
#'
#' @export
arima_predict <-
  function(model,test) {

    tm_ids <- c(1, grep("Tm", colnames(test)))
    test_xreg <- subset(test, select = -tm_ids)
    if (ncol(test_xreg) < 1) {
      test_xreg <- NULL
    } else {
      test_xreg <- as.matrix(test_xreg)
    }

    y <- get_y(test, model$form)

    preds <- as.vector(forecast::Arima(y = y,
                                       model=model,
                                       xreg = test_xreg)$fitted)
    # preds <- tryCatch(as.vector(forecast::Arima(y = y, model=model, xreg = test_xreg)$fitted),
    #                   error = function(e) NA)

    # if (is.na(preds[1])) {
    #   preds <- rep(mean(y), times = length(y))
    # }

    preds
  }


#' Fitting an Exponential Smoothing model
#' Using function \strong{ets} from forecast package
#'
#' @param y time series struture after applying
#' function \strong{embed_timeseries}
#'
#' @param form formula
#'
#' @keywords internal
#'
#' @export
bm_ets <-
  function(y, form) {
    y <- get_y(y, form)
    utils::capture.output(model <- forecast::ets(y = y))

    model$form <- form

    model
  }




#' Predict function for a model created from bm_ets
#'
#' @param model model created using \strong{bm_ets} function
#'
#' @param test test data created using \strong{embed_timeseries} function
#'
#' @keywords internal
#'
#' @export
ets_predict <-
  function(model,test) {
    y <- get_y(test, model$form)

    suppressMessages(preds <- tryCatch(forecast::ets(y = y, model = model),
                      error = function(e) NA))

    if (is.na(preds[1])) {
      trainset <- model$x
      predsf <- vnapply(1:length(y), function(j) {
        trainset <- c(trainset, y[seq_len(j) - 1])
        pforecast <- forecast::ets(trainset,1)$mean

        as.vector(pforecast)
      })
    } else {
      predsf <- as.vector(preds$fitted)
    }

    as.vector(predsf)
  }

#' Fitting an tbats model
#' Using function \strong{tbats} from forecast package
#'
#' @param y time series struture after applying
#' function \strong{embed_timeseries}
#'
#' @param form formula
#'
#' @keywords internal
#'
#' @export
bm_tbats <-
  function(y, form) {
    y <- get_y(y, form)
    utils::capture.output(model <- forecast::tbats(y = y))
    model$form <- form

    model
  }

#' Predict function for a model created from bm_tbats
#'
#' @param model model created using \strong{bm_tbats} function
#'
#' @param test test data created using \strong{embed_timeseries} function
#'
#' @keywords internal
#'
#' @export
tbats_predict <-
  function(model,test) {
    y <- get_y(test, model$form)

    preds <- tryCatch(as.vector(forecast::tbats(y = y, model = model)$fitted.values),
                      error = function(e) {
                        rep(mean(model$y, times = length(y)))
                      })

    preds
  }

#' Classical time series models
#'
#' Fit classical time series models, such as ARIMA
#' or exponential smoothing, to a time series
#'
#' @param form formula
#' @param data training data (see example)
#' @param lpars list of parameters (see example)
#'
#' @examples
#'
#' data("water_consumption")
#'
#' x <- embed_timeseries(water_consumption, 5)[1:100,]
#'
#' form <- target ~.
#' lpars <-
#' list(bm_timeseries = list(model = c("bm_arima","bm_tbats")))
#'
#' bmodels <- bm_timeseries(form, x, lpars)
#'
#' @keywords internal
#'
#' @export
bm_timeseries <-
  function(form, data, lpars) {
    Y <- get_y(data, form)

    if (is.null(lpars)) {
      lpars$bm_timeseries$model <- "bm_arima"
    }

    nmodels <- length(lpars$timeseries$model)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (modeltype in lpars$bm_timeseries$model) {
        j <- j + 1L

        mnames[j] <- gsub("bm_","",modeltype)

        utils::capture.output(ensemble[[j]] <-
          do.call(modeltype, list(y = data, form = form)))

    }
    names(ensemble) <- mnames

    ensemble
  }


#' Fit Gaussian Process models
#'
#' Learning a Gaussian Process model from training
#' data. Parameter setting can vary in \strong{kernel}
#' and \strong{tolerance}. See \code{\link[kernlab]{gausspr}}
#' for a comprehensive description.
#'
#' Imports learning procedure from \strong{kernlab} package.
#'
#' @param form formula
#' @param data training data for building the predictive
#' model
#' @param lpars a list containing the learning parameters
#'
#' @family base learning models
#'
#' @seealso other learning models: \code{\link{bm_mars}};
#' \code{\link{bm_ppr}}; \code{\link{bm_gbm}};
#' \code{\link{bm_glm}}; \code{\link{bm_cubist}};
#' \code{\link{bm_randomforest}}; \code{\link{bm_pls_pcr}};
#' \code{\link{bm_ffnn}}; \code{\link{bm_svr}}
#'
#' @import kernlab
#'
#' @return A list containing Gaussian Processes models
#'
#' @keywords internal
#'
#' @export
bm_gaussianprocess <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_gaussianprocess))
      lpars$bm_gaussianprocess <- list()

    if (is.null(lpars$bm_gaussianprocess$kernel))
      lpars$bm_gaussianprocess$kernel <- "rbfdot"
    if (is.null(lpars$bm_gaussianprocess$tol))
      lpars$bm_gaussianprocess$tol <- 0.001

    nmodels <-
      length(lpars$bm_gaussianprocess$kernel) *
      length(lpars$bm_gaussianprocess$tol)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (kernel in lpars$bm_gaussianprocess$kernel) {
      for (tolerance in lpars$bm_gaussianprocess$tol) {
        j <- j + 1L

        mnames[j] <- paste("gp", kernel, "kernel", tolerance, "tl", sep = "_")
        cat(mnames[j],"\n")
        if (!is.null(lpars$rm_ids)) {
          if (mnames[j] %in% names(lpars$rm_ids)) {
            rm_ids <- lpars$rm_ids[[mnames[j]]]
            data <- data[-rm_ids, ]
          }
        }

        ensemble[[j]] <-
          gausspr(form,
                  data,
                  type = "regression",
                  kernel = kernel,
                  tol = tolerance)

      }
    }
    names(ensemble) <- mnames

    ensemble
  }

#' Fit Projection Pursuit Regression models
#'
#' Learning a Projection Pursuit Regression
#' model from training data. Parameter setting
#' can vary in \strong{nterms} and \strong{sm.method}
#' parameters. See \code{\link[stats]{ppr}} for a comprehensive description.
#'
#' Imports learning procedure from \strong{stats} package.
#'
#' @inheritParams bm_gaussianprocess
#'
#' @family base learning models
#'
#' @seealso other learning models: \code{\link{bm_mars}};
#' \code{\link{bm_gaussianprocess}}; \code{\link{bm_gbm}};
#' \code{\link{bm_glm}}; \code{\link{bm_cubist}};
#' \code{\link{bm_randomforest}}; \code{\link{bm_pls_pcr}};
#' \code{\link{bm_ffnn}}; \code{\link{bm_svr}}
#'
#' @importFrom stats ppr
#' @keywords internal
#'
#' @export
bm_ppr <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_ppr))
      lpars$bm_ppr <- list()

    if (is.null(lpars$bm_ppr$nterms))
      lpars$bm_ppr$nterms <- 3
    if (is.null(lpars$bm_ppr$sm.method))
      lpars$bm_ppr$sm.method <- "supsmu"

    nmodels <-
      length(lpars$bm_ppr$nterms) * length(lpars$bm_ppr$sm.method)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (nterm in lpars$bm_ppr$nterms) {
      for (smoother in lpars$bm_ppr$sm.method) {
        j <- j + 1L
        mnames[j] <- paste0("ppr_", nterm, "_nterms_", smoother,"_method")
        cat(mnames[j],"\n")
        if (!is.null(lpars$rm_ids)) {
          if (mnames[j] %in% names(lpars$rm_ids)) {
            rm_ids <- lpars$rm_ids[[mnames[j]]]
            data <- data[-rm_ids, ]
          }
        }

        ensemble[[j]] <-
          ppr(form,
              data,
              nterms = nterm,
              sm.method = smoother)

      }
    }
    names(ensemble) <- mnames

    ensemble
  }

#' Fit Generalized Linear Models
#'
#' Learning a Generalized Linear Model
#' from training data. Parameter setting
#' can vary in \strong{alpha}.
#' See \code{\link[glmnet]{glmnet}} for a comprehensive description.
#'
#' Imports learning procedure from \strong{glmnet} package.
#'
#' @inheritParams bm_gaussianprocess
#'
#' @family base learning models
#'
#' @seealso other learning models: \code{\link{bm_mars}};
#' \code{\link{bm_ppr}}; \code{\link{bm_gbm}};
#' \code{\link{bm_gaussianprocess}}; \code{\link{bm_cubist}};
#' \code{\link{bm_randomforest}}; \code{\link{bm_pls_pcr}};
#' \code{\link{bm_ffnn}}; \code{\link{bm_svr}}
#'
#' @import glmnet
#' @keywords internal
#'
#' @export
bm_glm <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_glm))
      lpars$bm_glm <- list()

    if (is.null(lpars$bm_glm$alpha)) {
      lpars$bm_glm$alpha <- c(0, 1)
    }

    if (is.null(lpars$bm_glm$family)) {
      lpars$bm_glm$family <- "gaussian"
    }

    X <- model.matrix.na(form, data)
    Y <- get_y(data, form)

    nmodels <-
      length(lpars$bm_glm$alpha) *
      length(lpars$bm_glm$family)


    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (alpha in lpars$bm_glm$alpha) {
      for (fam in lpars$bm_glm$family) {
        j <- j + 1L

        if (alpha == 0) {
          mnames[j] <- "glm_ridge"
        } else if (alpha == 1) {
          mnames[j] <- "glm_lasso"
        } else {
          mnames[j] <- paste("glm_enet", alpha, sep = "_")
        }
        cat(mnames[j],"\n")

        mnames[j] <- paste(mnames[j],fam,sep="_")

        m.all <- glmnet(X, Y, alpha = alpha, family = fam)
        ensemble[[j]] <-
          glmnet(X,
                 Y,
                 alpha = alpha,
                 lambda = min(m.all$lambda),
                 family = fam)
      }
    }
    names(ensemble) <- mnames

    ensemble
  }


#' Fit Generalized Boosted Regression models
#'
#' Learning a Boosted Tree Model
#' from training data. Parameter setting
#' can vary in \strong{interaction.depth},
#' \strong{n.trees}, and \strong{shrinkage}
#' parameters.
#'
#' See \code{\link[gbm]{gbm}} for a comprehensive description.
#'
#' Imports learning procedure from \strong{gbm} package.
#'
#' @family base learning models
#'
#' @seealso other learning models: \code{\link{bm_mars}};
#' \code{\link{bm_ppr}}; \code{\link{bm_gaussianprocess}};
#' \code{\link{bm_glm}}; \code{\link{bm_cubist}};
#' \code{\link{bm_randomforest}}; \code{\link{bm_pls_pcr}};
#' \code{\link{bm_ffnn}}; \code{\link{bm_svr}}
#'
#' @inheritParams bm_gaussianprocess
#'
#' @import gbm
#' @keywords internal
#'
#' @export
bm_gbm <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_gbm))
      lpars$bm_gbm <- list()

    if (is.null(lpars$bm_gbm$interaction.depth))
      lpars$bm_gbm$interaction.depth <- 1
    if (is.null(lpars$bm_gbm$shrinkage))
      lpars$bm_gbm$shrinkage <- 0.001
    if (is.null(lpars$bm_gbm$n.trees))
      lpars$bm_gbm$n.trees <- 500
    if (is.null(lpars$bm_gbm$dist))
      lpars$bm_gbm$dist <- "gaussian"

    gbm_p <- lpars$bm_gbm
    nmodels <-
      length(gbm_p$interaction.depth) *
      length(gbm_p$shrinkage) *
      length(gbm_p$n.trees) *
      length(gbm_p$dist)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (id in lpars$bm_gbm$interaction.depth) {
      for (mdist in lpars$bm_gbm$dist) {
        for (shrinkage in lpars$bm_gbm$shrinkage) {
          for (n.trees in lpars$bm_gbm$n.trees) {
            j <- j + 1L

            mnames[j] <-
              paste("gbm", mdist, n.trees, "t", id, "id", shrinkage, "sh", sep = "_")

            cat(mnames[j],"\n")
            if (!is.null(lpars$rm_ids)) {
              if (mnames[j] %in% names(lpars$rm_ids)) {
                rm_ids <- lpars$rm_ids[[mnames[j]]]
                data <- data[-rm_ids, ]
              }
            }

            ensemble[[j]] <-
              gbm(
                form,
                data,
                distribution = mdist,
                interaction.depth = id,
                shrinkage = shrinkage,
                n.trees = n.trees
              )

          }
        }
      }
    }
    names(ensemble) <- mnames

    ensemble
  }

#' Fit Random Forest models
#'
#' Learning a Random Forest Model
#' from training data. Parameter setting
#' can vary in \strong{num.trees} and \strong{mtry}
#' parameters.
#'
#' See \code{\link[ranger]{ranger}} for a comprehensive description.
#'
#' Imports learning procedure from \strong{ranger} package.
#'
#' @family base learning models
#'
#' @seealso other learning models: \code{\link{bm_mars}};
#' \code{\link{bm_ppr}}; \code{\link{bm_gbm}};
#' \code{\link{bm_glm}}; \code{\link{bm_cubist}};
#' \code{\link{bm_gaussianprocess}}; \code{\link{bm_pls_pcr}};
#' \code{\link{bm_ffnn}}; \code{\link{bm_svr}}
#'
#' @inheritParams bm_gaussianprocess
#'
#' @import ranger
#' @keywords internal
#'
#' @export
bm_randomforest <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_randomforest))
      lpars$bm_randomforest <- list()

    if (is.null(lpars$bm_randomforest$num.trees))
      lpars$bm_randomforest$num.trees <- 500
    if (is.null(lpars$bm_randomforest$mtry))
      lpars$bm_randomforest$mtry <- ceiling(ncol(data) / 3)

    bad_mtry <- lpars$bm_randomforest$mtry > (ncol(data) - 1)

    if (any(bad_mtry)) {
      b_id <- which(bad_mtry)
      lpars$bm_randomforest$mtry[b_id] <- ceiling(ncol(data) / 3)
    }

    nmodels <-
      length(lpars$bm_randomforest$num.trees) *
      length(lpars$bm_randomforest$mtry)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (num.trees in lpars$bm_randomforest$num.trees) {
      for (mtry in lpars$bm_randomforest$mtry) {
        j <- j + 1L

        mnames[j] <- paste0("rf_n_", num.trees, "_m_", mtry)
        cat(mnames[j],"\n")
        if (!is.null(lpars$rm_ids)) {
          if (mnames[j] %in% names(lpars$rm_ids)) {
            rm_ids <- lpars$rm_ids[[mnames[j]]]
            data <- data[-rm_ids, ]
          }
        }


        ensemble[[j]] <-
          ranger(
            form,
            data,
            num.trees = num.trees,
            mtry = mtry,
            write.forest = TRUE)
      }
    }
    names(ensemble) <- mnames

    ensemble
  }

#' Fit Cubist models (M5)
#'
#' Learning a M5 model from training data
#' Parameter setting can vary in \strong{committees}
#' and \strong{neighbors} parameters.
#'
#' See \code{\link[Cubist]{cubist}} for a comprehensive description.
#'
#' Imports learning procedure from \strong{Cubist} package.
#'
#' @inheritParams bm_gaussianprocess
#'
#' @family base learning models
#'
#' @seealso other learning models: \code{\link{bm_mars}};
#' \code{\link{bm_ppr}}; \code{\link{bm_gbm}};
#' \code{\link{bm_glm}}; \code{\link{bm_gaussianprocess}};
#' \code{\link{bm_randomforest}}; \code{\link{bm_pls_pcr}};
#' \code{\link{bm_ffnn}}; \code{\link{bm_svr}}
#'
#' @importFrom Cubist cubist
#' @keywords internal
#'
#' @export
bm_cubist <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_cubist))
      lpars$bm_cubist <- list()

    if (is.null(lpars$bm_cubist$committees))
      lpars$bm_cubist$committees <- 50
    if (is.null(lpars$bm_cubist$neighbors))
      lpars$bm_cubist$neighbors <- 0

    form <- stats::as.formula(paste(deparse(form), "-1"))

    nmodels <-
      length(lpars$bm_cubist$committees) *
      length(lpars$bm_cubist$neighbors)

    X <- model.matrix.na(form, data)
    Y <- get_y(data, form)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (ncom in lpars$bm_cubist$committees) {
      for (neighbors in lpars$bm_cubist$neighbors) {
        j <- j + 1L

        mnames[j] <- paste0("cub_", ncom, "it", neighbors, "_nn")
        cat(mnames[j],"\n")

        if (!is.null(lpars$rm_ids)) {
          if (mnames[j] %in% names(lpars$rm_ids)) {
            rm_ids <- lpars$rm_ids[[mnames[j]]]
            X <- X[-rm_ids, ]
            Y <- Y[-rm_ids]
          }
        }


        ensemble[[j]] <-
          cubist(X, Y, committees = ncom, neighbors = neighbors)

        ensemble[[j]]$neighbors <- neighbors

      }
    }
    names(ensemble) <- mnames

    ensemble
  }

#' Fit Multivariate Adaptive Regression Splines models
#'
#' Learning a Multivariate Adaptive Regression Splines
#' model from training data.
#'
#' Parameter setting can vary in \strong{nk},
#' \strong{degree}, and \strong{thresh} parameters.
#'
#' See \code{\link[earth]{earth}} for a comprehensive description.
#'
#' Imports learning procedure from \strong{earth} package.
#'
#' @inheritParams bm_gaussianprocess
#'
#' @family base learning models
#'
#' @seealso other learning models: \code{\link{bm_gaussianprocess}};
#' \code{\link{bm_ppr}}; \code{\link{bm_gbm}};
#' \code{\link{bm_glm}}; \code{\link{bm_cubist}};
#' \code{\link{bm_randomforest}}; \code{\link{bm_pls_pcr}};
#' \code{\link{bm_ffnn}}; \code{\link{bm_svr}}
#'
#' @importFrom earth earth
#' @keywords internal
#'
#' @export
bm_mars <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_mars))
      lpars$bm_mars <- list()

    if (is.null(lpars$bm_mars$nk))
      lpars$bm_mars$nk <- 10
    if (is.null(lpars$bm_mars$degree))
      lpars$bm_mars$degree <- 3
    if (is.null(lpars$bm_mars$thresh))
      lpars$bm_mars$thresh <- 0.001

    if (is.null(lpars$bm_mars$pmethod))
      lpars$bm_mars$pmethod <- "backward"

    nmodels <-
      length(lpars$bm_mars$nk) *
      length(lpars$bm_mars$degree) *
      length(lpars$bm_mars$thresh) *
      length(lpars$bm_mars$pmethod)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (nk in lpars$bm_mars$nk) {
      for (method_ in lpars$bm_mars$pmethod) {
      for (degree in lpars$bm_mars$degree) {
        for (thresh in lpars$bm_mars$thresh) {
          j <- j + 1L

          mnames[j] <- paste0("mars_nk_", nk, "_d_", degree, "_t_", thresh,"_m_",method_)
          cat(mnames[j],"\n")
          if (!is.null(lpars$rm_ids)) {
            if (mnames[j] %in% names(lpars$rm_ids)) {
              rm_ids <- lpars$rm_ids[[mnames[j]]]
              data <- data[-rm_ids, ]
            }
          }


          ensemble[[j]] <-
            earth(form,
                  data,
                  nk = nk,
                  degree = degree,
                  thresh = thresh,
                  pmethod=method_)

        }
      }
      }
    }
    names(ensemble) <- mnames

    ensemble
  }


#' Fit Support Vector Regression models
#'
#' Learning a Support Vector Regression
#' model from training data.
#'
#' Parameter setting can vary in \strong{kernel},
#' \strong{C}, and \strong{epsilon} parameters.
#'
#' See \code{\link[kernlab]{ksvm}} for a comprehensive description.
#'
#' Imports learning procedure from \strong{kernlab} package.
#'
#' @inheritParams bm_gaussianprocess
#'
#' @family base learning models
#'
#' @seealso other learning models: \code{\link{bm_mars}};
#' \code{\link{bm_ppr}}; \code{\link{bm_gbm}};
#' \code{\link{bm_glm}}; \code{\link{bm_cubist}};
#' \code{\link{bm_randomforest}}; \code{\link{bm_pls_pcr}};
#' \code{\link{bm_ffnn}}; \code{\link{bm_gaussianprocess}}
#'
#' @import kernlab
#' @keywords internal
#'
#' @export
bm_svr <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_svr))
      lpars$bm_svr <- list()

    if (is.null(lpars$bm_svr$scale))
      lpars$bm_svr$scale <- FALSE
    if (is.null(lpars$bm_svr$type))
      lpars$bm_svr$type <- "eps-svr"
    if (is.null(lpars$bm_svr$kernel))
      lpars$bm_svr$kernel <- "vanilladot"
    if (is.null(lpars$bm_svr$epsilon))
      lpars$bm_svr$epsilon <- 0.1
    if (is.null(lpars$bm_svr$C))
      lpars$bm_svr$C <- 1

    nmodels <-
      length(lpars$bm_svr$kernel) *
      length(lpars$bm_svr$epsilon) *
      length(lpars$bm_svr$C)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (kernel in lpars$bm_svr$kernel) {
      for (epsilon in lpars$bm_svr$epsilon) {
        for (C in lpars$bm_svr$C) {
          j <- j + 1L

          mnames[j] <- paste0("svm_", kernel, "_g_", epsilon, "c_", C)
          cat(mnames[j],"\n")
          if (!is.null(lpars$rm_ids)) {
            if (mnames[j] %in% names(lpars$rm_ids)) {
              rm_ids <- lpars$rm_ids[[mnames[j]]]
              data <- data[-rm_ids, ]
            }
          }


          ensemble[[j]] <-
            ksvm(
              form,
              data,
              scale = lpars$bm_svr$scale,
              kernel = kernel,
              type = lpars$bm_svr$type,
              epsilon = epsilon,
              C = C)

        }
      }
    }
    names(ensemble) <- mnames

    ensemble
  }


#' Fit Feedforward Neural Networks models
#'
#' Learning a Feedforward Neural Network
#' model from training data.
#'
#' Parameter setting can vary in \strong{size}, \strong{maxit},
#' and \strong{decay} parameters.
#'
#' See \code{\link[nnet]{nnet}} for a comprehensive description.
#'
#' Imports learning procedure from \strong{nnet} package.
#'
#' @inheritParams bm_gaussianprocess
#'
#' @family base learning models
#'
#' @seealso other learning models: \code{\link{bm_mars}};
#' \code{\link{bm_ppr}}; \code{\link{bm_gbm}};
#' \code{\link{bm_glm}}; \code{\link{bm_cubist}};
#' \code{\link{bm_randomforest}}; \code{\link{bm_pls_pcr}};
#' \code{\link{bm_gaussianprocess}}; \code{\link{bm_svr}}
#'
#' @import monmlp
#' @keywords internal
#'
#' @export
bm_ffnn <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_ffnn))
      lpars$bm_ffnn <- list()

    if (is.null(lpars$bm_ffnn$hidden1))
      lpars$bm_ffnn$hidden1 <- 10
    if (is.null(lpars$bm_ffnn$hidden2))
      lpars$bm_ffnn$hidden2 <- 0

    nmodels <-
      length(lpars$bm_ffnn$hidden1) *
      length(lpars$bm_ffnn$hidden2)

    X <- as.matrix(model.matrix.na(form, data))
    Y <- as.matrix(get_y(data, form))

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (n1 in lpars$bm_ffnn$hidden1) {
        for (n2 in lpars$bm_ffnn$hidden2) {
          j <- j + 1L

          mnames[j] <- paste0("nnet_s1_", n1, "_s2_", n2)
          cat(mnames[j],"\n")

          ensemble[[j]] <-
            monmlp.fit(
              X,
              Y,
              hidden1=n1,
              hidden2=n2,
              n.ensemble=1,
              bag=F,
              silent=T
            )

      }
    }
    names(ensemble) <- mnames

    ensemble
  }


# bm_ffnn <-
#   function(form, data, lpars) {
#     if (is.null(lpars$bm_ffnn))
#       lpars$bm_ffnn <- list()
#
#     if (is.null(lpars$bm_ffnn$trace))
#       lpars$bm_ffnn$trace <- FALSE
#     if (is.null(lpars$bm_ffnn$linout))
#       lpars$bm_ffnn$linout <- TRUE
#     if (is.null(lpars$bm_ffnn$size))
#       lpars$bm_ffnn$size <- 30
#     if (is.null(lpars$bm_ffnn$decay))
#       lpars$bm_ffnn$decay <- 0.01
#     if (is.null(lpars$bm_ffnn$maxit))
#       lpars$bm_ffnn$maxit <- 750
#
#     nmodels <-
#       length(lpars$bm_ffnn$maxit) *
#       length(lpars$bm_ffnn$size) *
#       length(lpars$bm_ffnn$decay)
#
#     j <- 0
#     ensemble <- vector("list", nmodels)
#     mnames <- character(nmodels)
#     for (maxit in lpars$bm_ffnn$maxit) {
#       for (size in lpars$bm_ffnn$size) {
#         for (decay in lpars$bm_ffnn$decay) {
#           j <- j + 1L
#
#           mnames[j] <- paste0("nnet_s_", size, "_d_", decay, "_m_", maxit)
#           cat(mnames[j],"\n")
#           if (!is.null(lpars$rm_ids)) {
#             if (mnames[j] %in% names(lpars$rm_ids)) {
#               rm_ids <- lpars$rm_ids[[mnames[j]]]
#               data <- data[-rm_ids, ]
#             }
#           }
#
#           ensemble[[j]] <-
#             nnet(
#               form,
#               data,
#               linout = lpars$bm_ffnn$linout,
#               size = size,
#               maxit = maxit,
#               decay = decay,
#               trace = lpars$bm_ffnn$trace,
#               MaxNWts = 1000000
#             )
#
#         }
#       }
#     }
#     names(ensemble) <- mnames
#
#     ensemble
#   }


#' Fit PLS/PCR regression models
#'
#' Learning aPartial Least Squares or
#' Principal Components Regression from training data
#'
#' Parameter setting can vary in \strong{method}
#'
#' See \code{\link[pls]{mvr}} for a comprehensive description.
#'
#' Imports learning procedure from \strong{pls} package.
#'
#' @param form formula
#' @param data data to train the model
#' @param lpars parameter setting: For this multivariate regression
#' model the main parameter is "method". The available options are
#' "kernelpls", "svdpc", "cppls", "widekernelpls", and "simpls"
#'
#' @importFrom pls mvr
#'
#' @family base learning models
#'
#' @seealso other learning models: \code{\link{bm_mars}};
#' \code{\link{bm_ppr}}; \code{\link{bm_gbm}};
#' \code{\link{bm_glm}}; \code{\link{bm_cubist}};
#' \code{\link{bm_randomforest}}; \code{\link{bm_gaussianprocess}};
#' \code{\link{bm_ffnn}}; \code{\link{bm_svr}}
#'
#' @keywords internal
#'
#' @export
bm_pls_pcr <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_pls_pcr))
      lpars$bm_pls_pcr <- list()

    if (is.null(lpars$bm_pls_pcr$method))
      lpars$bm_pls_pcr$method <- "kernelpls"

    nmodels <-
      length(lpars$bm_pls_pcr$method)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (method in lpars$bm_pls_pcr$method) {
      j <- j + 1L

      mnames[j] <- paste("mvr", method, sep = "_")
      cat(mnames[j],"\n")

      if (!is.null(lpars$rm_ids)) {
        if (mnames[j] %in% names(lpars$rm_ids)) {
          rm_ids <- lpars$rm_ids[[mnames[j]]]
          data <- data[-rm_ids, ]
        }
      }

      model <-
        tryCatch(mvr(formula = form,
                     data = data,
                     method = method), error = function(e) NULL)

      if (!is.null(model)) {
        model$best_comp_train <- best_mvr(model, form, data)
      } else {
        mnames[j] <- NA_character_
      }
      ensemble[[j]] <- model
    }
    mnames <- mnames[!is.na(mnames)]
    ## se nalgum ensemble existe modelo, mete nomes...
    all_null <- all(sapply(ensemble, is.null))
    if (!all_null) { ### ensemble
      names(ensemble) <- mnames
    }

    ensemble
  }

#' Get best PLS/PCR model
#'
#' @param obj PLS/PCR model object
#' @param form formula
#' @param validation_data validation data used for
#' predicting performances of the model by number
#' of principal components
#'
#' @keywords internal
#'
#' @export
best_mvr <-
  function(obj, form, validation_data) {
    val_hat <- predict(obj, validation_data)

    target_var <- get_target(form)
    Y <- get_y(validation_data, form)

    val_hat <- as.data.frame(val_hat)

    err_by_comp <-
      sapply(val_hat,
             function(o)
               rmse(Y, o),
             USE.NAMES = FALSE)

    which.min(err_by_comp)
  }

#' predict method for pls/pcr
#'
#' @param model pls/pcr model
#' @param newdata new data
#'
#' @keywords internal
predict_pls_pcr <-
  function(model, newdata) {
    bcomp <- model$best_comp_train
    as.data.frame(predict(model, newdata))[,bcomp]
  }
