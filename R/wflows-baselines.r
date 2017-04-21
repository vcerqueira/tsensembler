#' stacking
#'
#' @inheritParams boostedADE
#'
#' @import ranger
#'
#' @export
Stacking <- function(form, train, test, learner, learner.pars, ...) {
  K <- ncol(train)
  target <- get_target(form)

  M <- learnM(form, train, learner, learner.pars, K)

  Y <- get_y(test, form)

  Y_hat <- predict(M, test)
  Y_hat.prop <- prop_hat(Y_hat)
  Y_hat_set <- cbind.data.frame(target = Y, Y_hat.prop)

  OOB.train <- ForwardValidation(x = train, nfolds = 10, OOB.fun.hat, .rbind = TRUE,
                                 form = form,
                                 learner = learner,
                                 learner.pars = learner.pars,
                                 embedding.dimension = K)

  oM <- ranger(target ~., OOB.train, num.trees = 1000, write.forest = TRUE)

  y_hat <- predict(oM, Y_hat_set)$predictions

  list(preds = y_hat, trues = Y)
}

#' Forecaster Combination
#' Model according to the guidelines
#' of Sanchez in adaptive combination of forecasts
#'
#' @param form Formula
#' @param train embedded time series used for training the base learners
#' @param test embedded time series used for testing
#' @param learner Character vector describing the base algorithms to be trained.
#' @param learner.pars Named list describing the parameter of the \code{learner}. Below are
#' described some examples.
#' @param ff forgetting factor
#' @param ... Further parameters to pass to the function
#'
#' @export
AEC <- function(form, train, test,
                learner, learner.pars, ff = .95, ...) {
  K <- get_embedsize(train)
  target <- get_target(form)
  Y <- unname(get_y(test, form))
  seq. <- seq_along(Y)

  M <- learnM(form, train, learner, learner.pars, K)
  
  Y_hat <- predict(M, test)
  Y_hat.prop <- prop_hat(Y_hat)
  e_y <- loss_M(Y_hat, lossFUN = ae)
  
  W <- as.data.frame(matrix(nrow = NROW(test),
                            ncol = ncol(Y_hat.prop)
                            )
                     )
  W[1:2, ] <- 1.
  
  Wy <- as.data.frame(lapply(seq_along(W), function(y_W) {
    y_hat <- Y_hat.prop[,y_W]
    E_y <- e_y[,y_W]

    r <- vnapply(seq_along(y_hat)[-(1:2)], function(y_w) {
      var_y <- var(y_hat[seq_len(y_w - 1)])
      v <- (1 / sqrt(var_y))

      v * exp(-((E_y[y_w - 1]^2) / (2*var_y))) * ff
    })
    c(1.,1., r)
  }))
  colnames(Wy) <- colnames(Y_hat.prop)
  
  bad_m <- which(sapply(Wy, function(j) any(is.na(j))))
  
  if (length(bad_m) > 0) {
    Wy <- Wy[,-bad_m]
    Y_hat.prop <- Y_hat.prop[,-bad_m]
  }

  y_hat <- vnapply(seq., function(j) {
    sum(unlistn(Y_hat.prop[j, ]) * proportion(unlistn(Wy[j, ])))
  })

  res <- list(trues = Y, preds = y_hat)

  res
}

#' Forecaster Combination sw50
#'
#' @param form Formula
#' @param train embedded time series used for training the base learners
#' @param test embedded time series used for testing
#' @param learner Character vector describing the base algorithms to be trained.
#' @param learner.pars Named list describing the parameter of the \code{learner}. Below are
#' described some examples.
#' @param ma_n moving avg
#' @param ... Further parameters to pass to the function
#'
#' @export
S_W50 <- function(form, train, test,
                  learner, learner.pars, ma_n = 50,  ...) {
  K <- get_embedsize(train)
  target <- get_target(form)
  Y <- unname(get_y(test, form))
  seq. <- seq_along(Y)

  M <- learnM(form, train, learner, learner.pars, K)

  Y_hat <- predict(M, test)
  Y_hat.prop <- prop_hat(Y_hat)

  E_Y <- loss_M(Y_hat, lossFUN = se)

  MASE <- rollmeanmatrix(E_Y, ma.N=ma_n)
  preweights <- rep(1./ncol(Y_hat.prop), times = ncol(Y_hat.prop))

  fScores <- data.frame(t(apply(MASE, 1, function(j) {
    proportion(normalizeMaxMin(-j))
  })))
  W <- rbind.data.frame(preweights, fScores[-NROW(fScores), ])

  y_hat <- vnapply(seq., function(j) {
    sum(unlistn(Y_hat.prop[j, ]) * proportion(unlistn(W[j, ])))
  })

  res <- list(trues = Y, preds = y_hat)

  res
}

#' Forecaster Combination
#' Model according to the guidelines
#' of Timmermann in Elusive Return Predictability
#'
#' @param form Formula
#' @param train embedded time series used for training the base learners
#' @param test embedded time series used for testing
#' @param learner Character vector describing the base algorithms to be trained.
#' @param learner.pars Named list describing the parameter of the \code{learner}. Below are
#' described some examples.
#' @param ma_n moving avg
#' @param R2_min min r squared
#' @param ... Further parameters to pass to the function
#'
#' @export
ERP_Timmermann <- function(form, train, test,
                           learner, learner.pars,
                           ma_n = 50, R2_min = .1, ...) {
  K <- get_embedsize(train)
  target <- get_target(form)
  Y <- unname(get_y(test, form))
  seq. <- seq_along(Y)

  M <- learnM(form, train, learner, learner.pars, K)

  Y_hat <- predict(M, test)
  Y_hat.prop <- prop_hat(Y_hat)

  R2 <- rbind_(lapply(seq., function(j) {
    if (j <= ma_n)
      in_seq <- seq_len(j)
    else
      in_seq <- (j-ma_n-1):j

    as.data.frame(
      lapply(Y_hat.prop[in_seq, ],
             function(y_hat)
               r_squared(Y[in_seq], y_hat)
      )
    )
  }))
  R2[1, ] <- 1

  condition1 <- apply(R2, 1, function(j) any(j > R2_min))

  prevailing_mean <- vnapply(seq.[-1], function(j) {
    if (j <= ma_n + 1)
      in_seq <- seq_len(j) - 1 #golden rule
    else
      in_seq <- (j-ma_n-1):j - 1 #golden rule

    mean(Y[in_seq])
  })

  prevailing_mean <- c(0, prevailing_mean) # no primeiro uso sempre combinacao

  y_hat <- vnapply(seq., function(j) {
    if (condition1[[j]])
      y_hat <- mean(unlistn(Y_hat.prop[j, ]))
    else
      y_hat <- prevailing_mean[j]

    y_hat
  })

  res <- list(trues = Y, preds = y_hat)

  res
}

#' blending
#'
#' @inheritParams boostedADE
#'
#' @import ranger
#'
#' @export
Blending <- function(form, train, test, learner, learner.pars, ...) {
  K <- ncol(train)
  target <- get_target(form)

  M <- learnM(form, train, learner, learner.pars, K)

  Y <- get_y(test, form)

  Y_hat <- predict(M, test)
  Y_hat.prop <- prop_hat(Y_hat)
  Y_hat_set <- cbind.data.frame(target = Y, Y_hat.prop)

  OOB.train <- holdout(x = train, estimation.ratio = .7, FUN = OOB.fun.hat,
                       form = form,
                       learner = learner,
                       learner.pars = learner.pars,
                       embedding.dimension = K)

  oM <- ranger(target ~., OOB.train, num.trees = 1000, write.forest = TRUE)
  y_hat <- predict(oM, Y_hat_set)$predictions

  list(preds = y_hat, trues = Y)
}

#' ss_W
#'
#' @param form form
#' @param train train
#' @param test test
#' @param learner learner
#' @param learner.pars learner pars
#' @param varying.embed ve
#' @param varying.trainwindow vw
#' @param ... adsa
#'
#' @export
s_w <- function(form, train, test, learner, learner.pars,
                    varying.embed = FALSE,
                    varying.trainwindow = FALSE, ...) {
  tgt <- get_target(form)
  Y <- get_y(test, form)

  embed.cols <- which(get_embedcols(train))

  embedding.dimension <- get_embedsize(train)
  row.size   <- NROW(train)

  #traintarget <- train[ ,tgt]
  #testtarget  <- test[ ,tgt]

  #train <- cbind.data.frame(target = traintarget, embedStats(train[ ,embed.cols]))
  #test  <- cbind.data.frame(target = testtarget,  embedStats(test[ ,embed.cols]))

  nstats <- ncol(train) - embedding.dimension
  if (varying.embed) {
    Ksplits <- c(embedding.dimension, embedding.dimension/2, embedding.dimension/3)
  } else {
    Ksplits <- embedding.dimension
  }
  if (varying.trainwindow) {
    Rsplits <- floor(c(row.size, row.size/2, row.size/4)) - 1
  } else {
    Rsplits <- row.size
  }

  Learningensemble <- Learntse(train,
                             form,
                             learner=learner,
                             learner.pars=learner.pars,
                             Ksplits = Ksplits,
                             Rsplits = Rsplits,
                             embedding.dimension=embedding.dimension,
                             nstats=0,
                             verbose=FALSE)

  .tse <- tseModel(baseModels = Learningensemble$learningModels,
                 preWeights = Learningensemble$preweights,
                 form = form,
                 ma.N = 30,
                 embedding.dimension = embedding.dimension,
                 committee.ratio = .1,
                 aggregationFUN = "static-s",
                 Mstar = Learningensemble$Mstar)

  Y_hat <- predict(.tse, test)
  Y_hat.prop <- prop_hat(Y_hat)

  Y_hat.tr <- predict(.tse, train)
  Y_hat.tr.prop <- prop_hat(Y_hat.tr)
  Y.tr <- get_y(train, form)

  n.test <- nrow(test)
  seq. <- seq_len(n.test)

  W <- list()
  ext.Y_hat <- Y_hat.tr.prop
  ext.Y <- Y.tr
  for (i in seq.) {
    W[[i]] <- apply(ext.Y_hat, 2, function(j) rmse(j, ext.Y))

    ext.Y_hat <- rbind.data.frame(ext.Y_hat, Y_hat.prop[seq_len(i), ])
    ext.Y <- c(ext.Y, Y[seq_len(i)])
  }

  W <- do.call(rbind.data.frame, W)
  colnames(W) <- colnames(Y_hat.tr.prop)

  W <- t(apply(W, 1, modelWeighting))

  y_hat <- vnapply(seq., function(j) sum(Y_hat.prop[j,] * W[j, ]))

  list(trues = Y, preds = y_hat)
}


#' ss_80
#'
#' @param form form
#' @param train train
#' @param test test
#' @param learner learner
#' @param learner.pars learner pars
#' @param varying.embed ve
#' @param varying.trainwindow vw
#' @param ... adsa
#'
#' @export
s_s80 <- function(form, train, test, learner, learner.pars,
                    varying.embed = TRUE,
                    varying.trainwindow = TRUE, ...) {
  tgt <- get_target(form)
  Y <- get_y(test, form)

  embed.cols <- which(get_embedcols(train))

  embedding.dimension <- length(embed.cols) + 1
  row.size   <- NROW(train)

  traintarget <- train[ ,tgt]
  testtarget  <- test[ ,tgt]

  train <- cbind.data.frame(target = traintarget, embedStats(train[ ,embed.cols]))
  test  <- cbind.data.frame(target = testtarget,  embedStats(test[ ,embed.cols]))

  nstats <- ncol(train) - embedding.dimension
  if (varying.embed) {
    Ksplits <- c(embedding.dimension, embedding.dimension/2, embedding.dimension/3)
  } else {
    Ksplits <- embedding.dimension
  }
  if (varying.trainwindow) {
    Rsplits <- floor(c(row.size, row.size/2, row.size/4)) - 1
  } else {
    Rsplits <- row.size
  }

  Learningensemble <- Learntse(train,
                             form,
                             learner = learner,
                             learner.pars = learner.pars,
                             Ksplits = Ksplits,
                             Rsplits = Rsplits,
                             embedding.dimension=embedding.dimension,
                             nstats=2,
                             verbose=FALSE)

  .tse <- tseModel(baseModels = Learningensemble$learningModels,
                 preWeights = Learningensemble$preweights,
                 form = form,
                 ma.N = 30,
                 embedding.dimension = embedding.dimension,
                 committee.ratio = .2,
                 aggregationFUN = "static-s",
                 Mstar = Learningensemble$Mstar)

  Y_hat <- predict(.tse, test)
  Y_hat.prop <- prop_hat(Y_hat)

  Y_hat.tr <- predict(.tse, train)
  Y_hat.tr.prop <- prop_hat(Y_hat.tr)
  Y.tr <- get_y(train, form)

  n.test <- nrow(test)
  seq. <- seq_len(n.test)

  W <- list()
  ext.Y_hat <- Y_hat.tr.prop
  ext.Y <- Y.tr
  for (i in seq.) {
    W[[i]] <- apply(ext.Y_hat, 2, function(j) rmse(j, ext.Y))

    ext.Y_hat <- rbind.data.frame(ext.Y_hat, Y_hat.prop[seq_len(i), ])
    ext.Y <- c(ext.Y, Y[seq_len(i)])
  }

  W <- do.call(rbind.data.frame, W)
  colnames(W) <- colnames(Y_hat.tr.prop)

  W <- t(apply(W, 1, modelWeighting))

  y_hat <- vnapply(seq., function(j) {
    w <- W[j, ]
    trimmedW <- which(w > quantile(w, .2))
    mean(unlist(Y_hat.prop[j, trimmedW]))
  })

  list(trues = Y, preds = y_hat)
}

#' ss_champ
#'
#' @param form form
#' @param train train
#' @param test test
#' @param learner learner
#' @param learner.pars learner pars
#' @param varying.embed ve
#' @param varying.trainwindow vw
#' @param ... adsad
#'
#' @export
s_champ <- function(form, train, test, learner, learner.pars,
                    varying.embed = TRUE,
                    varying.trainwindow = TRUE, ...) {
  tgt <- get_target(form)
  Y <- get_y(test, form)

  embed.cols <- which(get_embedcols(train))

  embedding.dimension <- length(embed.cols) + 1
  row.size   <- NROW(train)

  traintarget <- train[ ,tgt]
  testtarget  <- test[ ,tgt]

  train <- cbind.data.frame(target = traintarget, embedStats(train[ ,embed.cols]))
  test  <- cbind.data.frame(target = testtarget,  embedStats(test[ ,embed.cols]))

  nstats <- ncol(train) - embedding.dimension
  if (varying.embed) {
    Ksplits <- c(embedding.dimension, embedding.dimension/2, embedding.dimension/3)
  } else {
    Ksplits <- embedding.dimension
  }
  if (varying.trainwindow) {
    Rsplits <- floor(c(row.size, row.size/2, row.size/4)) - 1
  } else {
    Rsplits <- row.size
  }

  Learningensemble <- Learntse(train,
                               form,
                               learner=learner,
                               learner.pars=learner.pars,
                               Ksplits = Ksplits,
                               Rsplits = Rsplits,
                               embedding.dimension=embedding.dimension,
                               nstats=2,
                               verbose=FALSE)

  .tse <- tseModel(baseModels = Learningensemble$learningModels,
                   preWeights = Learningensemble$preweights,
                   form = form,
                   ma.N = 30,
                   embedding.dimension = embedding.dimension,
                   committee.ratio = .2,
                   aggregationFUN = "static-s",
                   Mstar = Learningensemble$Mstar)

  Y_hat <- predict(.tse, test)
  Y_hat.prop <- prop_hat(Y_hat)

  AE <- loss_M(Y_hat, prop = TRUE, ae)
  Y_hat.tr <- predict(.tse, train)
  Y_hat.tr.prop <- prop_hat(Y_hat.tr)
  Y.tr <- get_y(train, form)

  AE.tr <- loss_M(Y_hat.tr, prop = TRUE, ae)

  n.test <- nrow(test)
  seq. <- seq_len(n.test)

  W <- list()
  ext.AE <- AE.tr
  for (i in seq.) {
    w <- select_best(ext.AE)
    W[[i]] <- proportion(normalizeMaxMin(colSums(w)))

    ext.AE <- rbind.data.frame(ext.AE, AE[seq_len(i), ])
  }
  W <- do.call(rbind.data.frame, W)
  colnames(W) <- colnames(Y_hat.tr.prop)

  y_hat <- vnapply(seq., function(j) sum(Y_hat.prop[j,] * W[j, ]))

  list(trues = Y, preds = y_hat)
}

#' ss_kmeans
#'
#' @param form form
#' @param train train
#' @param test test
#' @param learner learner
#' @param learner.pars learner pars
#' @param varying.embed ve
#' @param varying.trainwindow vw
#' @param ... sadad
#'
#' @export
s_kmeans <- function(form, train, test, learner, learner.pars,
                     varying.embed = TRUE,
                     varying.trainwindow = TRUE, ...) {
  tgt <- get_target(form)
  Y <- get_y(test, form)

  embed.cols <- which(get_embedcols(train))

  embedding.dimension <- length(embed.cols) + 1
  row.size   <- NROW(train)

  traintarget <- train[ ,tgt]
  testtarget  <- test[ ,tgt]

  train <- cbind.data.frame(target = traintarget, embedStats(train[ ,embed.cols]))
  test  <- cbind.data.frame(target = testtarget,  embedStats(test[ ,embed.cols]))

  nstats <- ncol(train) - embedding.dimension
  if (varying.embed) {
    Ksplits <- c(embedding.dimension, embedding.dimension/2, embedding.dimension/3)
  } else {
    Ksplits <- embedding.dimension
  }
  if (varying.trainwindow) {
    Rsplits <- floor(c(row.size, row.size/2, row.size/4)) - 1
  } else {
    Rsplits <- row.size
  }

  Learningensemble <- Learntse(train,
                               form,
                               learner=learner,
                               learner.pars=learner.pars,
                               Ksplits = Ksplits,
                               Rsplits = Rsplits,
                               embedding.dimension=embedding.dimension,
                               nstats=2,
                               verbose=FALSE)

  .tse <- tseModel(baseModels = Learningensemble$learningModels,
                   preWeights = Learningensemble$preweights,
                   form = form,
                   ma.N = 30,
                   embedding.dimension = embedding.dimension,
                   committee.ratio = .2,
                   aggregationFUN = "static-s",
                   Mstar = Learningensemble$Mstar)

  Y_hat <- predict(.tse, test)
  Y_hat.prop <- prop_hat(Y_hat)

  SE <- loss_M(Y_hat, prop = TRUE, se)

  Y_hat.tr <- predict(.tse, train)
  Y_hat.tr.prop <- prop_hat(Y_hat.tr)
  Y.tr <- get_y(train, form)

  SE.tr <- loss_M(Y_hat.tr, prop = TRUE, se)

  n.test <- nrow(test)
  seq. <- seq_len(n.test)

  ext.SE <- SE.tr
  nK <- 3
  y_hat <- numeric(n.test)
  for (i in seq.) {
    kmeans.se <- kmeans(x = t(ext.SE), centers = 3, nstart = 10)

    ext.SE <- as.matrix(ext.SE)

    best_clust <- which.min(
      sapply(seq_len(nK), function(j) {
        clust.id <- which(kmeans.se$cluster == j)
        Reduce(mean, ext.SE[ , clust.id])
      })
    )

    bc.id <- which(kmeans.se$cluster == best_clust)

    y_hat[i] <- mean(unlist(Y_hat.prop[i, bc.id]))

    ext.SE <- rbind.data.frame(ext.SE, SE[seq_len(i), ])
  }

  list(trues = Y, preds = y_hat)
}
