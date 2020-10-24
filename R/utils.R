#' Get the response values from a data matrix
#'
#' Given a formula and a data set, \code{get_y} function retrieves
#' the response values.
#'
#' @param data data set with the response values;
#'
#' @param form formula
#'
#' @keywords internal
#'
#' @export
get_y <-
  function(data, form)
    stats::model.response(stats::model.frame(form, data, na.action = NULL))


#' Computing the rolling mean of the columns of a matrix
#'
#' @param x a numeric data.frame;
#' @param lambda periods to average over when computing the
#' moving average.
#'
#' @keywords internal
#'
#' @export
roll_mean_matrix <-
  function(x, lambda) {
    if (class(x) != "data.frame")
      stop("x obj class in roll_mean_matrix must be a data.frame")
    dim1 <- NROW(x)

    MASE <-
      lapply(x,
             function(z) {
               rollm <- RcppRoll::roll_mean(z, n = lambda)
               blanks <- dim1 - length(rollm)
               aux <- numeric(blanks)
               for (y in seq_along(aux)) {
                 aux[y] <- RcppRoll::roll_mean(z, n = y)[1]
               }
               c(aux, rollm)
             })
    as.data.frame(MASE)
  }

#' Splitting expressions by pattern
#'
#' This is an utility function that can be used to split expressions.
#' It is based on \code{strsplit} function.
#' split_by is the general purpose \emph{splitter}
#' split_by_ splits expressions by \"_\"
#' split_by. splits expressions by a dot
#'
#' @param expr character expression to split;
#'
#' @param split expression to split \code{expr} by;
#'
#' @param unlist. Logical. If TRUE, the splitted \code{expr} is unlisted;
#'
#' @param ... Further parameters to pass to \code{strsplit};
#'
#' @return a list or vector with a splitted expression
#'
#' @keywords internal
#'
#' @examples
#' split_by_("time_series")
#' split_by.("time.series")
#' split_by("born2bewild", "2")
#'
#' @export
split_by <- function(expr, split, unlist. = TRUE, ...) {
  expr <- strsplit(expr, split = split, fixed = TRUE, ...)
  if (unlist.) expr <- unlist(expr, use.names = FALSE)

  expr
}

#' @rdname split_by
#' @export
split_by_ <- function(expr, ...) split_by(expr, split = "_", ...)

#' @rdname split_by
#' @export
split_by. <- function(expr, ...) split_by(expr, split = ".", ...)


#' Scale a numeric vector using max-min
#'
#' Utility function used to linearly normalize a numeric vector
#'
#' @param x a numeric vector.
#'
#' @keywords internal
#'
#' @examples
#' normalize(rnorm(4L))
#' normalize(1:10)
#'
#' @return a linearly normalized vector
#'
#' @export
normalize <-
  function(x) {
    if (length(x) == 0L)
      stop("passed an argument of length 0.")
    if (!methods::is(x, "vector"))
      x <- unlist(x, use.names = FALSE)
    if (!is.numeric(x))
      stop("x must be numeric.")
    if (length(x) == 1L)
      return(1.)

    var_x <- stats::var(x, na.rm = T)
    if (is.na(var_x)) var_x <- 0

    if (var_x == 0) return(x)

    (x - min(x)) / (max(x) - min(x))
  }

#' Computing the proportions of a numeric vector
#'
#' Utility function used to compute the proportion of the
#' values of a vector. The proportion of a value is its ratio relative
#' to the sum of the vector.
#'
#' @param x a numeric vector;
#'
#' @examples
#' proportion(rnorm(5L))
#' proportion(1:10)
#'
#' @keywords internal
#'
#' @return A vector of proportions
#'
#' @export
proportion <-
  function(x) {
    if (!methods::is(x, "numeric"))
      x <- unlist(x)
    if (length(x) == 1L)
      return(1)

    var_x <- stats::var(x, na.rm = T)
    if (is.na(var_x)) var_x <- 0

    if (var_x == 0)
      return(rep(1 / length(x), times = length(x)))

    x / sum(x, na.rm = T)
  }

#' Get the target from a formula
#'
#' @param form formula
#'
#' @keywords internal
#'
#' @return the target variable as character
#'
#' @export
get_target <- function(form) split_by(deparse(form), " ")[1]

#' rbind with do.call syntax
#'
#' @param x object to call \code{rbind} to.
#'
#' @keywords internal
#'
#' @export
rbind_l <- function(x) do.call(rbind, x)

#' Applying lapply on the rows
#'
#' Wrapper function used to compute lapply on the rows
#' of a data.frame
#'
#' @param obj a data.frame object to apply the function.
#' @param FUN function to apply to each row of \code{obj}
#' @param ... Further parameters to \code{lapply}
#'
#' @keywords internal
#'
#' @export
l1apply <- function(obj, FUN, ...) {
  xseq <- seq_len(NROW(obj))
  lapply(X = xseq, FUN = FUN, ...)
}

#' Computing the softmax
#'
#' This function computes the softmax function in
#' a numeric vector
#'
#' @param x numeric vector
#'
#' @keywords internal
#'
#' @export
softmax <- function(x) {
  x <- x[!is.na(x)]

  exp(x) / sum(exp(x))
}

is_model_in_pars <-
  function(model, learner, lpnames) {
    par_list <- list(learner, lpnames)
    has_pars <-
      sapply(par_list,
             function(o)
               model %in% o)

    all(has_pars)
  }

are_pars_valid <-
  function(learner, lpars) {
    available_pars <-
      switch(learner,
             "bm_svr" = {
               c("kernel", "C", "epsilon")
             },
             "bm_glm" = {
               c("alpha","family")
             },
             "bm_gbm" = {
               c("interaction.depth", "shrinkage", "n.trees","dist")
             },
             "bm_ffnn" = {
               c("size", "decay", "maxit","hidden1","hidden2")
             },
             "bm_randomforest" = {
               c("num.trees", "mtry")
             },
             "bm_cubist" = {
               c("committees", "neighbors")
             },
             "bm_mars" = {
               c("nk", "thresh", "degree","pmethod")
             },
             "bm_timeseries" = {
               c("model")
             },
             "bm_gaussianprocess" = {
               c("kernel", "tol")
             },
             "bm_deepffnn" = {
               c("num_epochs", "nunits")
             },
             "bm_lstm" = {
               c("num_epochs", "nunits")
             },
             "bm_pls_pcr" = {
               c("method")
             },
             "bm_xgb" = {
               c("max_depth","eta","nrounds")
             },
             "bm_ppr" = {
               c("nterms", "sm.method")
             })

    mpars <- names(lpars[[learner]])

    if (!all(mpars %in% available_pars)) {
      bad_pars <- mpars[!mpars %in% available_pars]
      warning(paste("Bad or not supported specification of parameter",
                    "for model", learner, ":",
                    paste(bad_pars, collapse = ", ")))
    }
  }

extended_weights <-
  function(W, C) {
    stopifnot(NROW(W) == length(C))

    seq. <- seq_len(NROW(W))

    weightsf <-
      vapply(seq.,
             function(j) {
               W[j, C[[j]]] <- proportion(W[j, C[[j]]])
               W[j,-C[[j]]] <- 0.
               W[j, ]
             }, double(NCOL(W)))

    t(weightsf)
  }

subset_corr_matrix <-
  function(x, model) {
    corrs <- x[model,]
    corrs[!names(corrs) %in% model]
  }

get_embedcols <-
  function(x) grepl("^Tm[0-9]?[0-9]$", colnames(x))

rleid <-
  function(x) {
    x <- rle(x)$lengths
    rep(seq_along(x), times = x)
  }

rm.duplicated <-
  function(clist) {
    clist <-
      Map(function(x) paste0(x,collapse = "_"), clist)

    clist <- as.list(unique(unlist(clist)))

    Map(function(x) as.numeric(split_by_(x)), clist)
  }


model.matrix.na <-
  function(fm, x) {
    stats::model.matrix(fm,
                 stats::model.frame(fm,x,na.action=function(z) z))
  }
