#' Complementary Gaussian Error Function
#'
#' Erfc stands for the Complementary Gaussian Error Function.
#' This mathematical formula can be used as a squashing function.
#' Consider \code{x} a numeric vector representing the squared error of
#' base models in a given observation. By applying the erfc function on
#' the error, the weight of a given model decays exponentially as its
#' loss increases.
#'
#' @param x A numeric vector. The default value for the parameter
#' \code{lambda} presumes that \code{x} is in a 0--1 range. In the scope of
#' this package, this is achieved using normalize function;
#'
#' @param alpha parameter used to control the flatness of the erfc curve.
#' Defaults to 2.
#'
#' @return The complementary Gaussian error
#'
#' @references Cerqueira, Vitor; Torgo, Luis; Oliveira, Mariana,
#' and Bernhard Pfahringer. "Dynamic and Heterogeneous Ensembles
#' for Time Series Forecasting." Data Science and Advanced
#' Analytics (DSAA), 2017 IEEE International Conference on. IEEE, 2017.
#'
#' @examples
#' \dontrun{
#'   erfc(.1)
#'   erfc(c(.1, .7))
#' }
#'
#' @keywords internal
#'
#' @export
erfc <- function(x, alpha = 2) {
	stopifnot(is.numeric(x), is.numeric(alpha))

	2 * stats::pnorm(-sqrt(2) * x * alpha)
}

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

#' First-In First Out
#'
#' First-In First Out utility function inserts a
#' new value \code{inval} into a given sequential vector
#' \code{x}, dropping the last value of the sequence
#'
#' @param x a vector;
#'
#' @param inval new input value for vector x of the same
#' mode as \code{vector}
#'
#' @return A new vector \code{x}
#'
#' @keywords internal
#'
#' @examples
#' FIFO(1:10, 11)
#' FIFO(LETTERS[1:10], letters[1])
#'
#' @export
FIFO <-
  function(x, inval) {
    stopifnot(mode(x) == mode(inval))

    c(inval, x[-length(x)])
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
  if (unlist.) expr <- unlistn(expr)

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
#' @param ... Further arguments to min and max function
#' (e.g. na.rm = TRUE)
#'
#' @keywords internal
#'
#' @examples
#' normalize(rnorm(4L))
#' normalize(1:10)
#' normalize(c(1,2,NA,4), na.rm = TRUE)
#'
#' @return a linearly normalized vector
#'
#' @export
normalize <-
  function(x, ...) {
    if (length(x) == 0L)
      stop("passed an argument of length 0.")
    if (!methods::is(x, "vector"))
      x <- unlistn(x)
    if (!is.numeric(x))
      stop("x must be numeric.")
    if (length(x) == 1L)
      return(1.)
    if (stats::var(x, na.rm = T) == 0)
      return(x)

    (x - min(x, ...)) / (max(x, ...) - min(x, ...))
  }

#' Computing the proportions of a numeric vector
#'
#' Utility function used to compute the proportion of the
#' values of a vector. The proportion of a value is its ratio relative
#' to the sum of the vector.
#'
#' @param x a numeric vector;
#' @param ... Further arguments to \code{var} and \code{sum}
#' function (e.g. na.rm = TRUE)
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
  function(x, ...) {
    if (!methods::is(x, "numeric"))
      x <- unlist(x)
    if (length(x) == 1L)
      return(1)
    if (stats::var(x, na.rm = T) == 0)
      return(rep(1 / length(x), times = length(x)))

    x / sum(x, na.rm = T)
  }

#' vapply extension for logical values
#'
#' @param x A vector (atomic or list)
#' @param fun function to be applied
#' @param ... optional arguments to fun
#' @param use.names logical. Check argument USE.NAMES from \code{vapply} function.
#'
#' @seealso \code{\link{vapply}}
#'
#' @keywords internal
#'
#' @export
vlapply = function(x, fun, ..., use.names = FALSE) {
  vapply(X = x, FUN = fun, ..., FUN.VALUE = NA, USE.NAMES = use.names)
}

#' vapply extension for integer values
#'
#' @param x A vector (atomic or list)
#' @param fun function to be applied
#' @param ... optional arguments to fun
#' @param use.names logical. Check argument USE.NAMES from \code{vapply} function.
#'
#' @seealso \code{\link{vapply}}
#'
#' @keywords internal
#'
#' @export
viapply = function(x, fun, ..., use.names = FALSE) {
  vapply(X = x, FUN = fun, ..., FUN.VALUE = NA_integer_, USE.NAMES = use.names)
}

#' vapply extension for numeric values
#'
#' @param x A vector (atomic or list)
#' @param fun function to be applied
#' @param ... optional arguments to fun
#' @param use.names logical. Check argument USE.NAMES from \code{vapply} function.
#'
#' @seealso \code{\link{vapply}}
#'
#' @keywords internal
#'
#' @export
vnapply = function(x, fun, ..., use.names = FALSE) {
  vapply(X = x, FUN = fun, ..., FUN.VALUE = NA_real_, USE.NAMES = use.names)
}

#' vapply extension for character values
#'
#' @param x A vector (atomic or list)
#' @param fun function to be applied
#' @param ... optional arguments to fun
#' @param use.names logical. Check argument USE.NAMES from \code{vapply} function.
#'
#' @seealso \code{\link{vapply}}
#'
#' @keywords internal
#'
#' @export
vcapply = function(x, fun, ..., use.names = FALSE) {
  vapply(X = x, FUN = fun, ..., FUN.VALUE = NA_character_, USE.NAMES = use.names)
}

#' Unlist not using names
#'
#' @param x object to flat
#'
#' @seealso \code{unlist}
#'
#' @keywords internal
#'
#' @export
unlistn <- function(x) unlist(x, use.names = FALSE)

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

#' List without null elements
#'
#' @param l list with null elements
#'
#' @keywords internal
#'
#' @export
rm.null <- function(l) l[!vlapply(l, is.null)]


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
               c("size", "decay", "maxit")
             },
             "bm_randomforest" = {
               c("num.trees", "mtry")
             },
             "bm_cubist" = {
               c("committees", "neighbors")
             },
             "bm_mars" = {
               c("nk", "thresh", "degree")
             },
             "bm_gaussianprocess" = {
               c("kernel", "tol")
             },
             "bm_pls_pcr" = {
               c("method")
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

rm.len0 <-
  function(l) l[!vlapply(l, function(o) length(o) == 0)]

