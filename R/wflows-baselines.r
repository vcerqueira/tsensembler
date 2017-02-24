#' Dummy mean workflow
#' 
#' @param form Formula -- e.g. \code{target ~.}
#' @param train embedded time series for training
#' @param test embedded time series for testing
#' @param ... Further parameters to \code{meanhat} function
meanhat <- function(form, train, test, ...) {
  target <- get_target(form)
  preds <- mean(train[ ,target])

  list(preds = preds, trues = test[ ,target])
}

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
