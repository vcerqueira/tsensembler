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
        ensemble[[j]] <-
          gausspr(form,
                  data,
                  type = "regression",
                  kernel = kernel,
                  tol = tolerance)
        mnames[j] <- paste("gp", kernel, "krnl", sep = "_")
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
        ensemble[[j]] <-
          ppr(form,
              data,
              nterms = nterm,
              sm.method = smoother)

        mnames[j] <- paste0("ppr_", nterm, "nterms_", smoother)
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
#'
#' @export
bm_glm <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_glm))
      lpars$bm_glm <- list()

    if (is.null(lpars$bm_glm$alpha))
      lpars$bm_glm$alpha <- c(0, 1)

    X <- stats::model.matrix(form, data)
    Y <- get_y(data, form)

    nmodels <- length(lpars$bm_glm$alpha)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (alpha in lpars$bm_glm$alpha) {
      j <- j + 1L

      m.all <- glmnet(X, Y, alpha = alpha)
      ensemble[[j]] <-
        glmnet(X,
               Y,
               alpha = alpha,
               lambda = min(m.all$lambda))

      if (alpha == 0) {
        mnames[j] <- "glm_ridge"
      } else if (alpha == 1) {
        mnames[j] <- "glm_lasso"
      } else {
        mnames[j] <- paste("glm_enet", alpha, sep = "_")
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
      lpars$bm_gbm$n.trees <- 100

    gbm_p <- lpars$bm_gbm
    nmodels <-
      length(gbm_p$interaction.depth) *
      length(gbm_p$shrinkage) *
      length(gbm_p$n.trees)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (id in lpars$bm_gbm$interaction.depth) {
      for (shrinkage in lpars$bm_gbm$shrinkage) {
        for (n.trees in lpars$bm_gbm$n.trees) {
          j <- j + 1L

          ensemble[[j]] <-
            gbm(
              form,
              data,
              distribution = "gaussian",
              interaction.depth = id,
              shrinkage = shrinkage,
              n.trees = n.trees
            )

          mnames[j] <-
            paste("gbm", n.trees, "t", id, "id", shrinkage, "sh", sep = "_")
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
#'
#' @export
bm_randomforest <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_randomforest))
      lpars$bm_randomforest <- list()

    if (is.null(lpars$bm_randomforest$num.trees))
      lpars$bm_randomforest$num.trees <- 100
    if (is.null(lpars$bm_randomforest$mtry))
      lpars$bm_randomforest$mtry <- ncol(data) / 3

    nmodels <-
      length(lpars$bm_randomforest$num.trees) *
      length(lpars$bm_randomforest$mtry)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (num.trees in lpars$bm_randomforest$num.trees) {
      for (mtry in lpars$bm_randomforest$mtry) {
        j <- j + 1L
        ensemble[[j]] <-
          ranger(
            form,
            data,
            num.trees = num.trees,
            mtry = mtry,
            write.forest = TRUE)

        mnames[j] <- paste0("rf_n", num.trees, "m", mtry)
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
#'
#' @export
bm_cubist <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_cubist))
      lpars$bm_cubist <- list()

    if (is.null(lpars$bm_cubist$committees))
      lpars$bm_cubist$committees <- 30
    if (is.null(lpars$bm_cubist$neighbors))
      lpars$bm_cubist$neighbors <- 0

    form <- stats::as.formula(paste(deparse(form), "-1"))

    nmodels <-
      length(lpars$bm_cubist$committees) *
      length(lpars$bm_cubist$neighbors)

    X <- stats::model.matrix(form, data)
    Y <- get_y(data, form)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (ncom in lpars$bm_cubist$committees) {
      for (neighbors in lpars$bm_cubist$neighbors) {
        j <- j + 1L
        ensemble[[j]] <-
          cubist(X, Y, committees = ncom, neighbors = neighbors)

        mnames[j] <- paste0("cub_", ncom, "it", neighbors, "nn")
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
#'
#' @export
bm_mars <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_mars))
      lpars$bm_mars <- list()

    if (is.null(lpars$bm_mars$nk))
      lpars$bm_mars$nk <- 15
    if (is.null(lpars$bm_mars$degree))
      lpars$bm_mars$degree <- 3
    if (is.null(lpars$bm_mars$thresh))
      lpars$bm_mars$thresh <- 0.001

    nmodels <-
      length(lpars$bm_mars$nk) *
      length(lpars$bm_mars$degree) *
      length(lpars$bm_mars$thresh)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (nk in lpars$bm_mars$nk) {
      for (degree in lpars$bm_mars$degree) {
        for (thresh in lpars$bm_mars$thresh) {
          j <- j + 1L
          ensemble[[j]] <-
            earth(form,
                  data,
                  nk = nk,
                  degree = degree,
                  thresh = thresh)
          mnames[j] <- paste0("mars_nk", nk, "_d", degree, "t", thresh)
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

          ensemble[[j]] <-
            ksvm(
              form,
              data,
              scale = lpars$bm_svr$scale,
              kernel = kernel,
              type = lpars$bm_svr$type,
              epsilon = epsilon,
              C = C)

          mnames[j] <- paste0("svm_", kernel, "g", epsilon, "c", C)
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
#' @import nnet
#'
#' @export
bm_ffnn <-
  function(form, data, lpars) {
    if (is.null(lpars$bm_ffnn))
      lpars$bm_ffnn <- list()

    if (is.null(lpars$bm_ffnn$trace))
      lpars$bm_ffnn$trace <- FALSE
    if (is.null(lpars$bm_ffnn$linout))
      lpars$bm_ffnn$linout <- TRUE
    if (is.null(lpars$bm_ffnn$size))
      lpars$bm_ffnn$size <- 30
    if (is.null(lpars$bm_ffnn$decay))
      lpars$bm_ffnn$decay <- 0.01
    if (is.null(lpars$bm_ffnn$maxit))
      lpars$bm_ffnn$maxit <- 200

    nmodels <-
      length(lpars$bm_ffnn$maxit) *
      length(lpars$bm_ffnn$size) *
      length(lpars$bm_ffnn$decay)

    j <- 0
    ensemble <- vector("list", nmodels)
    mnames <- character(nmodels)
    for (maxit in lpars$bm_ffnn$maxit) {
      for (size in lpars$bm_ffnn$size) {
        for (decay in lpars$bm_ffnn$decay) {
          j <- j + 1L

          ensemble[[j]] <-
            nnet(
              form,
              data,
              linout = lpars$bm_ffnn$linout,
              size = size,
              maxit = maxit,
              decay = decay,
              trace = lpars$bm_ffnn$trace,
              MaxNWts = 1000000
            )

          mnames[j] <- paste0("nnet_s", size, "_d", decay, "m")
        }
      }
    }
    names(ensemble) <- mnames

    ensemble
  }


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

      ensemble[[j]] <-
        mvr(formula = form,
            data = data,
            method = method)

      ensemble[[j]]$best_comp_train <-
        best_mvr(ensemble[[j]], form, data)

      mnames[j] <- paste("mvr", method, sep = "_")
    }
    names(ensemble) <- mnames

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
