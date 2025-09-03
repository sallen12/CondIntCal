#' Interval Forecast Coverage
#'
#' @description
#' Calculate the (marginal) coverage of forecasts in the form of prediction intervals.
#'
#' @inheritParams is_decomp
#' @param closed logical specifying whether the coverage of the closed (\code{TRUE}) or open
#'  (\code{FALSE}) prediction intervals should be returned; default is \code{TRUE}
#' @param twosided logical specifying whether to return both the lower and upper
#'  exceedance probabilities rather than the combined coverage; default is \code{FALSE}
#'
#'
#' @return
#' The marginal coverage of the prediction intervals.
#'
#'
#' @details
#'
#' \emph{Theory:}
#'
#' Interval forecasts (or prediction intervals) are comprised of a lower bound \eqn{\ell}
#' and an upper bound \eqn{u}, with \eqn{\ell < u}. The forecast is made such that the
#' observation \eqn{y} is predicted to fall within the interval with a given coverage level \eqn{1 - \alpha}.
#' In the general case, it can be assumed that the prediction interval is \emph{non-central},
#' so that the probability that \eqn{y < \ell} is equal to \eqn{\alpha_1 \in (0, 1)} and the probability
#' that \eqn{y > u} is equal to \eqn{\alpha_2 \in (0, 1)}, with \eqn{\alpha_1 < \alpha_2}.
#' Typically, a \emph{central} prediction interval is issued, for which it is assumed that the
#' probability that \eqn{y < \ell} is equal to the probability that \eqn{y > u}, i.e.
#' \eqn{\alpha_1 = \alpha/2} and \eqn{\alpha_2 = 1 - \alpha/2}.
#'
#' In practice, we observe several interval forecasts \eqn{[\ell_i, u_i]} and corresponding observations \eqn{y_i}
#' for \eqn{i = 1, \dots, n}. The reliability of the interval forecasts can be assessed by checking whether the
#' empirical coverage is roughly equal to the nominal coverage level. That is, whether
#' \deqn{\frac{1}{n} \sum_{i=1}^{n} 1\{\ell_i \le y_i \le u_i\} \approxeq 1 - \alpha.}
#'
#' While this empirical coverage check is widely used in practice, it does not assess
#' whether the probability that \eqn{y < l} is equal to \eqn{\alpha_1 \in (0, 1)}, and the probability
#' that \eqn{y > u} is equal to \eqn{\alpha_2 \in (0, 1)}. This can be assessed by separately checking whether
#' \deqn{\frac{1}{n} \sum_{i=1}^{n} 1\{y_i < \ell_i\} \approxeq \alpha_1, \quad \frac{1}{n} \sum_{i=1}^{n} 1\{y_i > u_i\} \approxeq \alpha_2.}
#'
#'
#' \emph{Implementation:}
#'
#' \code{y} should be a numeric vector of observations.
#'
#' \code{int} should be a numeric matrix or dataframe containing the prediction intervals for \code{y}.
#' We should have that \code{nrow(int) = length(y)} and \code{ncol(int) = 2}. The first and
#' second columns of \code{int} should contain the lower and upper bounds of the prediction
#' intervals, respectively, so that \code{all(int[, 1] < int[, 2])}.
#'
#' \code{closed} should be a single logical. \code{closed = TRUE} (the default) corresponds to
#' \deqn{\frac{1}{n} \sum_{i=1}^{n} 1\{\ell_i \le y_i \le u_i\}}
#' while \code{closed = FALSE} returns
#' \deqn{\frac{1}{n} \sum_{i=1}^{n} 1\{\ell_i < y_i < u_i\}.}
#' The difference is only relevant in cases where there is a chance that the observation is
#' exactly equal to either the lower or upper bound of the interval forecast. This may be
#' the case when predicting a discrete or censored outcome variable, but otherwise the
#' two values will generally be the same.
#'
#' \code{twosided} should be a single logical. If \code{twosided = TRUE}, then rather
#' than returning the coverage, the individual exceedance probabilities of the lower and
#' upper interval bounds are returned. That is,
#' \deqn{\frac{1}{n} \sum_{i=1}^{n} 1\{y_i < \ell_i\}, \quad \frac{1}{n} \sum_{i=1}^{n} 1\{y_i > u_i\}.}
#' If \code{closed = FALSE}, then these are replaced with
#' \deqn{\frac{1}{n} \sum_{i=1}^{n} 1\{y_i \le \ell_i\}, \quad \frac{1}{n} \sum_{i=1}^{n} 1\{y_i \ge u_i\}.}
#' The coverage can be calculated using one minus the sum of the two returned exceedance probabilities.
#'
#'
#' @seealso \code{\link{is_decomp}} \code{\link{interval_score}}
#'
#'
#' @section References:
#'
#' Allen, S., Burnello, J. and Ziegel, J. (2025):
#' `Assessing the conditional calibration of interval forecasts using decompositions of the interval score'.
#' \emph{arXiv pre-print}.
#'
#'
#' @author Sam Allen
#'
#'
#' @examples
#' n <- 10000 # sample size
#' mu <- rnorm(n)
#' y <- rnorm(n, mean = mu, sd = 1) # simulate observations
#'
#' alpha <- 0.1 # 90% prediction intervals
#'
#' # Ideal forecaster: F = N(mu, 1)
#' L_id <- qnorm(alpha/2, mu)
#' U_id <- qnorm(1 - alpha/2, mu)
#' int_id <- data.frame(Lower = L_id, Upper = U_id)
#'
#' coverage(y, int_id)
#' coverage(y, int_id, twosided = TRUE)
#'
#'
#' @name coverage
NULL

# coverage
#' @export
#' @rdname coverage
coverage <- function(y, int, closed = TRUE, twosided = FALSE) {
  check_cov_args(y, int, closed, twosided)
  n <- length(y)
  if (n == 1) {
    if (!is.vector(int)) {
      y <- rep(y, nrow(int))
    } else {
      int <- int |> as.matrix() |> t()
    }
  } else {
    if (is.vector(int)) {
      int <- matrix(int, nrow = n, ncol = 2, byrow = T)
    }
  }
  if (closed) {
    ncov_low <- y < int[, 1]
    ncov_upp <- y > int[, 2]
  } else {
    ncov_low <- y <= int[, 1]
    ncov_upp <- y >= int[, 2]
  }
  if (twosided) {
    c("Lower" = mean(ncov_low), "Upper" = mean(ncov_upp))
  } else {
    1 - mean(ncov_low | ncov_upp)
  }
}

# argument checks
check_cov_args <- function(y, int, closed, twosided) {
  # y
  if (!is.numeric(y) || !is.vector(y)) stop("'y' must be a numeric vector")
  n <- length(y)

  # int
  if (!is.matrix(int) && !is.data.frame(int) && !is.vector(int)) stop("'int' must either be a vector, matrix or dataframe")
  if ((is.matrix(int) || is.vector(int)) && !is.numeric(int)) stop("'int' must contain numeric values")
  if (is.data.frame(int) && !all(sapply(int, is.numeric))) stop("'int' must contain numeric values")
  if (is.matrix(int) || is.data.frame(int)) {
    if (ncol(int) != 2) stop("'int' must be a matrix or dataframe of dimension (n, 2)")
    if (!(nrow(int) %in% c(1, n))) stop("the number of rows in 'int' must be either 1 or the same as the length of 'y'")
    if (any(int[, 2] < int[, 1])) stop("the second column of 'int' must always be no smaller than the first column")
  } else {
    if (length(int) != 2) stop("'int' must either be a vector of length 2, or a matrix or dataframe of dimension (n, 2)")
    if (int[2] < int[1]) stop("the second value of 'int' must be no smaller than the first")
  }

  # closed
  if (!is.logical(closed) || length(closed) > 1) stop("'closed' must be a single logical")

  # twosided
  if (!is.logical(twosided) || length(twosided) > 1) stop("'twosided' must be a single logical")
}
