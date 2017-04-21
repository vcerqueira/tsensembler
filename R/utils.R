#' Erfc
#'
#' Erfc stands for the Complementary Gaussian Error Function.
#' This mathematical formula can be used as a squashing function.
#'
#' @param x A numeric vector in the range 0-1
#' @param alpha parameter used to sharp the Erfc curve. Defaults to 1.
#'
#' @return The complementary Gaussian error
#'
#' @examples
#'
#' \dontrun{
#'   erfc(.1)
#'   erfc(c(.1, .7))
#' }
#'
#' @export
erfc <- function(x, alpha = 2) {
	stopifnot(is.numeric(x), is.numeric(alpha))

	2 * stats::pnorm(-sqrt(2) * x * alpha)
}

#' Auxiliary function used to compute column-wise moving averages on a matrix
#'
#' @param mat matrix
#' @param ma.N periods to average over when computing the moving average.
#'
#' @export
rollmeanmatrix <- function(mat, ma.N) {
  if (class(mat) != "data.frame") stop("mat obj class in rollmeanmatrix must be a data.frame")
  dim1 <- NROW(mat)
  MASE <- as.data.frame(
    lapply(mat, function(z) {
      rollm <- RcppRoll::roll_mean(z, n = ma.N)
      blanks <- dim1 - length(rollm)
      aux <- numeric(blanks)
      for (y in seq_along(aux)) {
        aux[y] <- RcppRoll::roll_mean(z, n = y)[1]
      }
    c(aux, rollm)
  }))
  MASE
}

#' First-In First Out
#'
#' First-In First Out utility function inserts a new value \code{.in} into
#' a given sequential vector \code{x}, dropping the last value of the sequence
#'
#' @param x Vector
#' @param .in New input value for vector x of the same class as \code{class(vector)}
#'
#' @return A new vector \code{x}
#'
#' @examples
#' FIFO(1:10, 11)
#' FIFO(LETTERS[1:10], letters[1])
#'
#' @export
FIFO <- function(x, .in) {
  stopifnot(mode(x) == mode(.in))

  c(.in, x[-length(x)])
}

#' Pull yourself up by your bootstraps
#'
#' bootstrap utility function for random sampling with replacement.
#'
#' @param n Vector to choose elements from and also the size of the sample
#'
#' @examples
#' \dontrun{
#' bootstrap("a")
#' bootstrap(FALSE)
#' }
#' bootstrap(5)
#' bootstrap(5)
#' bootstrap(3)
#'
#' @return A vector of sampled id's
#'
#' @export
bootstrap <- function(n) {
  if (!is.numeric(n)) stop("n must be numeric.")

  sample(n, n, replace = TRUE)
}

#' Splitting expressions by pattern
#'
#' This is an utility function that can be used to split expressions.
#' It is based on \code{strsplit} function.
#' split_by is the general purpose \emph{splitter}
#' split_by_ splits expressions by \"_\"
#' split_by. splits expressions by a dot
#'
#' @param expr character expression to split
#' @param split expression to split \code{expr} by
#' @param unlist. Logical. If TRUE, the splitted \code{expr} is unlisted.
#' @param ... Further parameters to pass to \code{strsplit}
#' @return a vector with a splitted expression
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

#' @rdname .splitBy
#' @export
split_by_ <- function(expr, ...) .splitBy(expr, split = "_", ...)

#' @rdname .splitBy
#' @export
split_by. <- function(expr, ...) .splitBy(expr, split = ".", ...)

#' normalize
#'
#' Utility function used to linearly normalize a numeric vector 
#' 
#' @param x a numeric vector.
#' @param ... Further arguments to min and max function (e.g. na.rm = TRUE)
#'
#' @examples
#' normalize(rnorm(4L))
#' normalize(1:10)
#' normalize(c(1,2,NA,4), na.rm = TRUE)
#'
#' @return a linearly normalized vector
#'
#' @export
normalize <- function(x, ...) {
  if (length(x) == 0L) {
    stop("passed an argument of length 0 on normalizeMaxMin function.")
  }
  if (length(x) == 1L) return(1.)
  if (!methods::is(x, "vector")) x <- unlistn(x)
  if (var(x, ...) == 0) return(x)
  if (!is.numeric(x)) stop("x must be numeric.")


  (x - min(x, ...)) / (max(x, ...) - min(x, ...))
}

#' proportion
#'
#' Utility function used to compute the proportion of the values of a vector
#' The proportion of a value is its ratio relative
#' to the sum of the vector.
#'
#' @param x a numeric vector.
#' @param ... Further arguments to \code{var} and \code{sum} function (e.g. na.rm = TRUE)
#'
#' @examples
#' proportion(rnorm(5L))
#' proportion(1:10)
#'
#' @return A vector of proportions
#'
#' @export
proportion <- function(x, ...) {
  if (!methods::is(x, "numeric")) x <- unlistn(x)
  if (length(x) == 1L) return(1)
  if (var(x, ...) == 0) return(rep(1/length(x), times = length(x)))

  x / sum(x, ...)
}


#' split_vec
#'
#' Utility function that splits a vector into n parts
#' In this package this function is used for smoothing large vectors for plotting
#'
#' @param x vector to split
#' @param n number of parts to split x into
#' @param avg Logical. If \code{TRUE} returns the results averaged by mean.
#'
#' @examples
#' split_vec(letters[1:6], 3, F)
#'
#' @return List of splitted vectors
#'
#' @export
split_vec <- function(x, n, avg = TRUE) {
  stopifnot(is.numeric(n))

  s_vec <- split(x, ceiling(seq_along(x) / n))

  if (avg) s_vec <- sapply(s_vec, mean)

  s_vec
}

#' Exponential Weighted Average Loss
#'
#' This is an utility function that computes the exponential
#' loss of the base learners, in order to measure the regret of the combined
#' model to such base models. This is based on theorem 2.3 presented by
#' CESA-BIANCHI and LUGOSI in their book Prediction, Learning, and Games
#' \(2006\).
#'
#' @param loss loss base models in a given prediction time.
#' @param N number of total experts
#' @param t_i prediction time. This is used in order to maintain uniform
#' loss bounds.
#' @param .proportion if TRUE, scales the results.
#' @param windsize Window size used to compute bounds
#'
#' @return Exponential loss of each base model in a given prediction time
#' step.
expLoss <- function(loss, N, t_i, .proportion = FALSE, windsize = NULL) {
  .p <- ifelse(is.null(windsize), t_i, windsize)
  zeta <- sqrt((8. * log(N)) / .p)

  .exp <- exp(-zeta * loss)
  if (.proportion) .exp <- proportion(.exp)
  .exp
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
#' @export
unlistn <- function(x) unlist(x, use.names = FALSE)

#' Get the target from a formula
#'
#' @param form formula
#'
#' @return the target variable as character
#'
#' @export
get_target <- function(form) .splitBy(deparse(form), " ")[1]


#' rbind with do.call syntax
#'
#' @param x object to call \code{rbind} to.
#'
#' @export
rbind_ <- function(x) do.call(rbind, x)


#' List without null elements
#'
#' @param l list with null elements
#'
#' @export
rm.null <- function(l) l[!vlapply(l, is.null)]

cmplt_rollmedian <- function(x, n) {
  len <- length(x)
  rollm <- RcppRoll::roll_median(x, n = n)
  blanks <- len - length(rollm)
  aux <- numeric(blanks)
  for (y in seq_along(aux)) {
    aux[y] <- RcppRoll::roll_median(x, n = y)[1]
  }
  c(aux, rollm)
}

#' rleid
#' 
#' @param x vector
rleid <- function(x) {
  x <- rle(x)$lengths
  rep(seq_along(x), times=x)
}

#' lapply on the rows
#'
#' @param obj obj
#' @param FUN func
#' @param ... Further parameters to \code{lapply}
#'
#' @export
l1apply <- function(obj, FUN, ...) {
  xseq <- seq_len(NROW(obj))
  lapply(X = xseq, FUN = FUN, ...)
}

#' diff log transformation of time series
#' 
#' @param x a timeseries
#' 
#' @export
diff_log <- function(x) {
  x <- diff(x)[-1]
  
  x.log <- sign(x) * log(abs(x) + 1)
  
  x.log
}

#' Incremental mean along a vector
#' 
#' @param x numeric vector
#' 
#' @export
incremental_mean <- function(x) {
  seq. <- seq_along(x)
  vnapply(seq., function(j) {
    sub_x <- x[seq_len(j)]
    RcppRoll::roll_mean(sub_x, length(sub_x))
  })
}

#' Incremental variance along a vector
#' 
#' @param x numeric vector
#' 
#' @export
incremental_var <- function(x) {
  seq. <- seq_along(x)
  vnapply(seq., function(j) {
    sub_x <- x[seq_len(j)]
    RcppRoll::roll_var(sub_x, length(sub_x))
  })
}

#' Logarithmic transformation
#' 
#' @param x numeric vector
#' 
#' @export
log_trans <- function(x) sign(x) * log(abs(x) + 1)

#' Softmax function
#' 
#' @param x numeric vector
#' 
#' @return squashed values according to the softmax
#' 
#' @export
softmax <- function(x) {
  x <- x[!is.na(x)]

  exp(x) / sum(exp(x))
}
