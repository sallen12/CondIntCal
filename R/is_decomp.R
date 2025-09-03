#' Interval Score Decompositions
#'
#' @description
#' Assess the conditional calibration of interval forecasts using decompositions of
#' the interval score.
#'
#' @param y vector of observations.
#' @param int matrix or dataframe containing the prediction intervals.
#' @param level nominal coverage level of the central prediction intervals.
#' @param alpha1,alpha2 lower and upper quantile levels if evaluating non-central prediction intervals.
#' @param method method used to estimate the decomposition terms; one of 'linear' (default) and 'isotonic'.
#' @param return_fit logical specifying whether the recalibrated predictions used to estimate
#'  the decomposition terms should be returned; default is \code{FALSE}.
#'
#'
#' @return
#'
#' Vector containing the average interval score for the interval forecasts,
#' as well as the estimated uncertainty, discrimination, and miscalibration terms.
#'
#' If \code{return_fit = TRUE}, a list is returned containing the interval score
#' decomposition terms and the recalibrated interval forecasts used to estimated
#' the decomposition terms.
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
#' Competing interval forecasts can be compared using the interval score,
#' \deqn{\mathrm{IS}_{\alpha_1, \alpha_2}([\ell, u], y) = |u - y | + \frac{1}{\alpha_1} 1\{y < \ell\} (\ell - y) + \frac{1}{1 - \alpha_2} 1\{y > u\} (y - u)}.
#' In the case of central prediction intervals, the scaling factors \eqn{1 / \alpha_1} and \eqn{1 / (1 - \alpha_2)} both
#' simplify to \eqn{2 / \alpha}.
#'
#' In practice, we observe several interval forecasts \eqn{[\ell_i, u_i]} and corresponding observations \eqn{y_i}
#' for \eqn{i = 1, \dots, n}, and we wish to compare forecasters or forecast methods based on
#' the average interval score
#' \deqn{\mathrm{\bar{S}} = \frac{1}{n} \sum_{i=1}^{n} \mathrm{IS}_{\alpha_1, \alpha_2}([\ell_i, u_i], y_i).}
#'
#' While the average score \eqn{\mathrm{\bar{S}}} provides a single value that can be used to rank
#' different forecasters, it can also be additively decomposed into terms quantifying different
#' aspects of forecast performance. Most commonly, these terms measure how much \emph{uncertainty} there
#' is in the forecast problem, the forecast's \emph{discrimination} or information content, and
#' the degree of forecast \emph{miscalibration}. These terms can be calculated using the following decomposition:
#' \deqn{\mathrm{\bar{S}} = \mathrm{\bar{S}^R} - (\mathrm{\bar{S}^R} - \mathrm{\bar{S}^C}) + (\mathrm{\bar{S}} - \mathrm{\bar{S}^C}),}
#' where
#' \deqn{\mathrm{\bar{S}^R} = \frac{1}{n} \sum_{i=1}^{n} \mathrm{IS}_{\alpha_1, \alpha_2}(R, y_i)}
#' is the average interval score for a reference prediction interval \eqn{R},
#' and
#' \deqn{\mathrm{\bar{S}^C} = \frac{1}{n} \sum_{i=1}^{n} \mathrm{IS}_{\alpha_1, \alpha_2}(C_i, y_i)}
#' is the average interval score for \eqn{C_i}, which correspond to recalibrated versions of
#' the original prediction intervals \eqn{[\ell_i, u_i]}.
#'
#' The first term of the decomposition, \eqn{\mathrm{\bar{S}^R}}, provides a baseline measure of forecast performance,
#' thereby quantifying the inherent uncertainty or unpredictability of the observations; the second term,
#' \eqn{\mathrm{\bar{S}^R} - \mathrm{\bar{S}^C}}, represents the improvement in accuracy that is obtained
#' from the recalibrated predictions compared to the (uninformative) reference predictions, thus providing a measure
#' of how much information is contained within the original interval forecasts; the third term,
#' \eqn{\mathrm{\bar{S}} - \mathrm{\bar{S}^C}}, is the improvement in accuracy that is obtained by recalibrating
#' the original interval forecasts, which quantifies the degree of miscalibration in the original predictions.
#' Note that the miscalibration corresponds to conditional calibration, where the conditioning is on the original
#' forecasts themselves.
#'
#' This decomposition holds for any choice of the reference prediction \eqn{R} and recalibration method
#' to produce \eqn{C_i}. However, the interpretation of the terms requires specific choices. In particular,
#' it is desirable that the decomposition terms are non-negative; in this case, we can say that a miscalibration
#' of zero corresponds to a calibrated prediction interval, for example.
#'
#' Two methods have been proposed in the literature to calculate the decomposition terms with
#' desirable theoretical guarantees. These correspond to recalibrating the original prediction intervals using
#' \emph{isotonic regression} and \emph{linear regression}. Details can be found in the references below.
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
#' \code{level} should be the nominal coverage of the prediction intervals. For example, if \code{int}
#' contains central 90% prediction intervals for \code{y}, then \code{level = 0.9}.
#' This corresponds to \eqn{1 - \alpha} in the theoretical explanation above.
#' If the prediction intervals are non-central prediction intervals for \code{y}, then
#' \code{alpha1} and \code{alpha2} should be used instead of \code{level}; if both are used,
#' then \code{level} will be ignored.
#'
#' \code{alpha1} and \code{alpha2} should be the lower and upper quantile levels corresponding
#' to the lower and upper bounds of \code{int}, if \code{int} contains non-central prediction
#' intervals for \code{y}. These correspond to \eqn{\alpha_1} and \eqn{\alpha_2} in the theoretical
#' explanation above. If the prediction intervals are central prediction intervals for \code{y}, then
#' it is simpler to use \code{level} instead of \code{alpha1} and \code{alpha2}.
#'
#' \code{method} should be the method used to recalibrate the original prediction intervals to estimate
#' the interval score decomposition terms. This must be one of \code{'isotonic'} and \code{'linear'},
#' corresponding to isotonic and linear regression, respectively.
#'
#' \code{return_fit} should be a logical specifying whether the recalibrated prediction intervals
#' should be returned alongside the decomposition terms. It may be useful to study the
#' recalibrated prediction intervals, in order to assess whether the assumptions made by
#' the recalibration method are appropriate, for example. The default is \code{return_fit = FALSE}.
#'
#'
#' @seealso \code{\link{interval_score}} \code{\link{coverage}} \code{\link{count_comparables}} \code{\link{plot_mcbdsc}}
#'
#'
#' @section References:
#'
#' \emph{Linear decomposition:}
#'
#' Dimitriadis, T., Puke, M., (2025).
#' `Statistical Inference for Score Decompositions under Linear Recalibration'
#'
#' \emph{Isotonic decomposition:}
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
#' n <- 1000 # sample size
#' alpha <- 0.1 # 90% prediction intervals
#' mu <- rnorm(n)
#' y <- rnorm(n, mean = mu, sd = 1) # simulate observations
#'
#' # Ideal forecaster: F = N(mu, 1)
#' L_id <- qnorm(alpha/2, mu)
#' U_id <- qnorm(1 - alpha/2, mu)
#' int_id <- data.frame(Lower = L_id, Upper = U_id)
#'
#' # isotonic decomposition
#' out_iso <- is_decomp(y, int_id, level = 1 - alpha) # using level
#' out_iso2 <- is_decomp(y, int_id, alpha1 = alpha/2, alpha2 = 1 - alpha/2) # using alpha1 and alpha2
#' out_iso3 <- is_decomp_iso(y, int_id, level = 1 - alpha) # using is_decomp_iso
#' all.equal(out_iso, out_iso2)
#' all.equal(out_iso, out_iso3)
#'
#' # linear decomposition
#' out_lin <- is_decomp(y, int_id, level = 1 - alpha, method = "linear") # using level
#' out_lin2 <- is_decomp(y, int_id, alpha1 = alpha/2, alpha2 = 1 - alpha/2, method = "linear") # using alpha1 and alpha2
#' out_lin3 <- is_decomp_lin(y, int_id, level = 1 - alpha) # using is_decomp_lin
#' all.equal(out_lin, out_lin2)
#' all.equal(out_lin, out_lin3)
#'
#' # non-central prediction intervals
#'
#' alpha1 <- 0.05
#' alpha2 <- 0.9
#'
#' # Ideal forecaster: F = N(mu, 1)
#' L_id <- qnorm(alpha1, mu)
#' U_id <- qnorm(alpha2, mu)
#' int_id <- data.frame(Lower = L_id, Upper = U_id)
#'
#' out_iso <- is_decomp(y, int_id, alpha1 = alpha1, alpha2 = alpha2) # isotonic
#' out_lin <- is_decomp(y, int_id, alpha1 = alpha1, alpha2 = alpha2, method = "linear") # linear
#'
#'
#' @name is_decomp
NULL

# interval score decomposition wrapper
#' @export
#' @rdname is_decomp
is_decomp <- function(y, int, level = NULL, alpha1 = NULL, alpha2 = NULL, method = c("isotonic", "linear"), return_fit = F) {
  method <- match.arg(method)

  if (method == "isotonic") {
    is_decomp_iso(y, int, level, alpha1, alpha2, return_fit)
  } else if (method == "linear") {
    is_decomp_lin(y, int, level, alpha1, alpha2, return_fit)
  }
}

# isotonic interval score decomposition
#' @export
#' @rdname is_decomp
is_decomp_iso <- function(y, int, level = NULL, alpha1 = NULL, alpha2 = NULL, return_fit = F) {
  check_args(y, int, level, alpha1, alpha2, method = "isotonic", return_fit)
  if (is.null(alpha1) && is.null(alpha2)) alpha1 <- (1 - level)/2; alpha2 <- 1 - alpha1

  int <- data.frame(Lower = int[, 1], Upper = int[, 2])

  idr_settings <- osqp::osqpSettings(verbose = F, eps_abs = 1e-10, eps_rel = 1e-5, max_iter=10000)
  fit <- isodistrreg::idr(y = y, X = int, pars = idr_settings) # fit IDR
  out <- stats::predict(fit) # get predicted distributions

  int_rc <- cbind(isodistrreg::qpred(out, alpha1), isodistrreg::qpred(out, alpha2)) # get recalibrated interval forecasts
  int_mg <- c(stats::quantile(y, alpha1, type = 1), stats::quantile(y, alpha2, type = 1)) # get unconditional interval forecasts

  IS <- interval_score(y, int, alpha1 = alpha1, alpha2 = alpha2) |> mean() # interval score of original forecast
  IS_mg <- interval_score(y, int_mg, alpha1 = alpha1, alpha2 = alpha2) |> mean() # interval score of unconditional forecast
  IS_rc <- interval_score(y, int_rc, alpha1 = alpha1, alpha2 = alpha2) |> mean() # interval score of recalibrated forecast

  MCB <- IS - IS_rc # miscalibration
  DSC <- IS_mg - IS_rc # discrimination

  if (isTRUE(return_fit)) {
    int_rc <- data.frame(Lower = int_rc[, 1], Upper = int_rc[, 2])
    return(list(decomp = c("IS" = IS, "UNC" = IS_mg, "DSC" = DSC, "MCB" = MCB), int_rc = int_rc))
  } else {
    return(c("IS" = IS, "UNC" = IS_mg, "DSC" = DSC, "MCB" = MCB))
  }
}

# linear interval score decomposition
#' @export
#' @rdname is_decomp
is_decomp_lin <- function(y, int, level = NULL, alpha1 = NULL, alpha2 = NULL, return_fit = F) {
  check_args(y, int, level, alpha1, alpha2, method = "linear", return_fit)
  if (is.null(alpha1) && is.null(alpha2)) alpha1 <- (1 - level)/2; alpha2 <- 1 - alpha1

  dat <- data.frame(Y = y, Lower = int[, 1], Upper = int[, 2])

  fit <- tryCatch(quantreg::rq(Y ~ Lower + Upper, data = dat, tau = c(alpha1, alpha2)),
                  error = function(e) quantreg::rq(Y ~ Lower, data = dat, tau = c(alpha1, alpha2)))
  int_rc <- predict(fit) # get recalibrated interval forecasts
  int_mg <- c(stats::quantile(y, alpha1, type = 1), stats::quantile(y, alpha2, type = 1)) # get unconditional interval forecasts

  IS <- interval_score(y, int, alpha1 = alpha1, alpha2 = alpha2) |> mean() # interval score of original forecast
  IS_mg <- interval_score(y, int_mg, alpha1 = alpha1, alpha2 = alpha2) |> mean() # interval score of unconditional forecast
  IS_rc <- interval_score(y, int_rc, alpha1 = alpha1, alpha2 = alpha2) |> mean() # interval score of recalibrated forecast

  MCB <- IS - IS_rc # miscalibration
  DSC <- IS_mg - IS_rc # discrimination

  if (isTRUE(return_fit)) {
    int_rc <- data.frame(Lower = int_rc[, 1], Upper = int_rc[, 2])
    return(list(decomp = c("IS" = IS, "UNC" = IS_mg, "DSC" = DSC, "MCB" = MCB), int_rc = int_rc))
  } else {
    return(c("IS" = IS, "UNC" = IS_mg, "DSC" = DSC, "MCB" = MCB))
  }
}

# argument checks
check_args <- function(y, int, level, alpha1, alpha2, method, return_fit) {
  # y
  if (!is.numeric(y) || !is.vector(y)) stop("'y' must be a numeric vector")
  n <- length(y)
  if (!(n > 1)) stop("'y' must contain multiple observations")

  # int
  if (!is.matrix(int) && !is.data.frame(int)) stop("'int' must either be a matrix or dataframe")
  if (is.matrix(int) && !is.numeric(int)) stop("'int' must contain numeric values")
  if (is.data.frame(int) && !all(sapply(int, is.numeric))) stop("'int' must contain numeric values")
  if (ncol(int) != 2) stop("'int' must be a matrix or dataframe of dimension (n, 2)")
  if (nrow(int) != n) stop("the number of rows in 'int' must be the same as the length of 'y'")
  if (any(int[, 2] < int[, 1])) stop("the second column of 'int' must always be no smaller than the first column")

  # level, alpha1, alpha2
  if (is.null(alpha1) && is.null(alpha2)) {
    # level
    if (is.null(level)) stop("at least one of 'level', 'alpha1', and 'alpha2' must be provided")
    if (!is.numeric(level) || length(level) > 1) stop("'level' must be a single numeric")
    if (level <= 0 || level >= 1) stop("'level' must be between 0 and 1")
  } else {
    if (!is.null(level)) warning("'level' and 'alpha1' and/or 'alpha2' are provided - 'level' will be ignored")
    # alpha 1
    if (!is.null(alpha1)) {
      if (!is.numeric(alpha1) || length(alpha1) > 1) stop("'alpha1' must be a single numeric")
      if (alpha1 <= 0 || alpha1 >= 1) stop("'alpha1' must be between 0 and 1")
    }
    # alpha 2
    if (!is.null(alpha2)) {
      if (!is.numeric(alpha2) || length(alpha2) > 1) stop("'alpha2' must be a single numeric")
      if (alpha2 <= 0 || alpha2 >= 1) stop("'alpha2' must be between 0 and 1")
    }
    if (alpha2 <= alpha1) stop("'alpha1' must be smaller than 'alpha2'")
  }

  # method
  if (!(method %in% c("isotonic", "linear"))) stop("'method' must be one of c('isotonic', 'linear')")
  if (method == "isotonic") {
    if (n < 500)
      warning(paste0("method = 'isotonic' typically requires at least 500 values to converge (here ", n, "). Consider using method = 'linear' instead"))
  }

  # return_fit
  if (!is.logical(return_fit) || length(return_fit) > 1) stop("'return_fit' must be a single logical")
}

