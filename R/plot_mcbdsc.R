#' Miscalibration-Discrimination Plot
#'
#' @description
#' Plot the miscalibration and discrimination ability of prediction methods, as
#' calculated from score decompositions. The average score, as well as the forecast
#' uncertainty, are also displayed.
#'
#' @param decomp A dataframe or list containing the decomposition terms.
#' @param n_isolines The number of isolines showing average scores.
#' @param colour_values A colour specification passed to the `values` argument
#'   of `scale_colour_manual()`. Recycled if length 1.
#' @param colour_unc A colour specification highlighting the UNC component layers.
#' @param MCBDSC_repel A logical specifying whether labels should be placed by the `ggrepel` package.
#' @param MCB_lim The plot limits for the x-axis (the MCB component).
#' @param DSC_lim The plot limits for the y-axis (the DSC component).
#'
#'
#' @return
#' The proportion of prediction intervals that are comparable (i.e. ordered).
#'
#'
#' @details
#'
#' The \code{plot_mcbdsc()} function is an adaptation of the \code{autoplot.triptych_mcbdsc()}
#' function in the `triptych` package.
#'
#' One way to visualise the output of score decompositions is to plot the
#' estimated miscalibration against the estimated discrimination. This can be done
#' for multiple prediction methods simultaneously. These so-called
#' \emph{miscalibration-discrimination plots} (or MCB-DSC plots) allow practitioners to compare
#' forecasters via their average score (using isolines on the plot),
#' whilst simultaneously analysing their calibration and information content.
#'
#' \code{decomp} should be a dataframe with \code{decomp$UNC} containing the uncertainty
#' terms of the score decompositions, \code{decomp$DSC} containing the discrimination terms,
#' \code{decomp$MCB} the miscalibration terms, and \code{decomp$forecast} containing labels
#' for the different prediction models.
#' Alternatively, \code{decomp} can be a list of named vectors with \code{UNC}, \code{DSC},
#' and \code{MCB} components. In this case, the list is automatically converted into a dataframe
#' of the above form, with the names of the list used as the labels for the different prediction
#' models. This list format works easily with the ouput of \code{\link{is_decomp}}, for example.
#'
#' \code{n_isolines} is the number of isolines shown on the plot, showing the average score
#' for different combinations of miscalibration and discrimination. This should be a single integer.
#'
#' \code{colour_values} should be a vector of colours, with length equal to the number of
#' prediction models (i.e. \code{nrow(decomp)}). If only one colour is provided, then
#' this is used for all methods.
#'
#' \code{colour_unc} should be a single colour. This corresponds to the colour used to
#' display the uncertainty term of the score decomposition.
#'
#' \code{MCBDSC_repel} should be a single logical. The default is \code{FALSE}, in which
#' case labels are not placed using the `ggrepel` package.
#'
#' \code{MCB_lim} and \code{DSC_lim} should be vectors of length two, containing the
#' lower and upper bounds of the x- and y-axes respectively. If only a single value is
#' given, then this is assumed to be the upper bound, with zero taken as the lower bound.
#'
#'
#' @seealso \code{\link{is_decomp}}
#'
#'
#' @section References:
#'
#' Dimitriadis, T., Gneiting, T., Jordan, A. I., & Vogel, P. (2024):
#' `Evaluating probabilistic classifiers: The triptych.'
#' \emph{International Journal of Forecasting}, 40, 1101-1122.
#'
#'
#' @examples
#' \dontrun{
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
#' out_iso <- is_decomp(y, int_id, level = 1 - alpha) # isotonic decomposition
#' out_lin <- is_decomp(y, int_id, level = 1 - alpha, method = "linear") # linear decomposition
#'
#' plot_mcbdsc(list(Isotonic = out_iso, Linear = out_lin), MCB_lim = out_iso[1], DSC_lim = out_iso[1])
#' }
#'
#' @name plot_mcbdsc
NULL

# MCB-DSC plot (adapted from triptych package)
#' @export
#' @rdname plot_mcbdsc
plot_mcbdsc <- function(decomp, n_isolines = 10, colour_values = "black", colour_unc = "#00BF7D", MCBDSC_repel = FALSE, MCB_lim = NA, DSC_lim = NA) {
  check_plot_args(decomp, n_isolines, colour_values, colour_unc, MCBDSC_repel, MCB_lim, DSC_lim)
  if (is.list(decomp)) {
    decomp <- decomp |> as.data.frame() |> t() |> as.data.frame()
    decomp$forecast <- rownames(decomp)
  }

  # Limits and out-of-bound identification
  default_lims <- function(x) c(0, 1.1 * max(x[is.finite(x)]))
  is_in_range <- function(x, xrange) x >= xrange[1] & x <= xrange[2]
  split_by_out_of_bounds <- function(df, MCB_lim, DSC_lim) {
    oob_state <- factor(
      x = with(df, dplyr::case_when(
        is_in_range(DSC, DSC_lim) & is_in_range(MCB, MCB_lim) ~ "within",
        is_in_range(DSC, DSC_lim) & !is.finite(MCB)           ~ "infty",
        TRUE                                                  ~ "oob"
      )),
      levels = c("within", "infty", "oob")
    )
    res <- split(df, oob_state)
    if (nrow(res$infty)) {
      res$infty$x_geom_text <- MCB_lim[2]
    }
    res
  }
  if (anyNA(MCB_lim)) {
    MCB_lim <- default_lims(decomp$MCB)
  } else if (length(MCB_lim) == 1) {
    MCB_lim <- c(0, MCB_lim)
  }
  if (anyNA(DSC_lim)) {
    DSC_lim <- default_lims(decomp$DSC)
  } else if (length(DSC_lim) == 1) {
    DSC_lim <- c(0, DSC_lim)
  }
  decomp_by_state <- split_by_out_of_bounds(decomp, MCB_lim, DSC_lim)
  # Check that the plot is not empty of points!
  if (!nrow(decomp_by_state$within)) {
    warning(paste(
      "The given limits for the MCB-DSC plot exclude all forecasts.",
      "The default choices are used instead."
    ))
    MCB_lim <- default_lims(decomp$MCB)
    DSC_lim <- default_lims(decomp$DSC)
    decomp_by_state <- split_by_out_of_bounds(decomp, MCB_lim, DSC_lim)
  }
  if (nrow(decomp_by_state$oob)) {
    message(paste(
      "The following forecasts are not included in the MCB-DSC plot as their",
      "miscalibration measure is outside the plot limits:",
      paste(decomp_by_state$oob$forecast, collapse = ", ")
    ))
  }

  # Reasonable score values for isolines
  choose_isolines <- function(unc, MCB_lim, DSC_lim) {
    scores <- pretty(x = unc - c(-1.1 * MCB_lim[2], DSC_lim[2]), n = n_isolines)
    tibble::tibble(
      slope = 1,
      intercept = unc - scores,
      label = scores
    ) |>
      # Remove a line if its intercept is too close to zero (less than 1/5 times the line distance).
      dplyr::filter(abs(.data$intercept) > abs(diff(.data$intercept)[1]) / 5)
  }
  df_iso_abline <- choose_isolines(decomp$UNC[1], MCB_lim, DSC_lim)

  # replicate colour_values if it has length 1
  if (length(colour_values) == 1L & nrow(decomp) > 1L) {
    colour_values <- rep(colour_values, nrow(decomp))
  }

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(x = 0, y = 0, xend = .data$max_val, yend = .data$max_val),
      data = tibble::tibble(max_val = 2 * max(MCB_lim, DSC_lim)),
      colour = colour_unc,
      linewidth = 1
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(.data$x, .data$y),
      data = data.frame(x = 0, y = 0),
      colour = colour_unc,
      fill = colour_unc,
      size = 2,
      shape = 15
    ) +
    ggplot2::geom_abline(
      data = df_iso_abline,
      mapping = ggplot2::aes(intercept = .data$intercept, slope = .data$slope),
      colour = "gray"
    ) +
    geomtextpath::geom_labelabline(
      data = df_iso_abline,
      mapping = ggplot2::aes(intercept = .data$intercept, slope = .data$slope, label = .data$label),
      colour = "gray",
      hjust = 0.85,
      size = 7 * 0.36,
      text_only = TRUE,
      boxcolour = NA,
      straight = TRUE
    ) +
    geomtextpath::geom_labelabline(
      mapping = ggplot2::aes(
        intercept = 0,
        slope = 1,
        label = paste("UNC:", prettyNum(.data$UNC, digits = 3))
      ),
      data = decomp[1, ],
      colour = colour_unc,
      hjust = 0.85,
      size = 7 * 0.36,
      text_only = TRUE,
      boxcolour = NA,
      straight = TRUE
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = .data$MCB, y = .data$DSC, colour = .data$forecast),
      data = decomp_by_state$within
    ) +
    {
      if (!isTRUE(MCBDSC_repel)) {
        ggplot2::geom_text(
          data = decomp_by_state$within,
          mapping = ggplot2::aes(x = .data$MCB, y = .data$DSC, label = .data$forecast, colour = .data$forecast),
          size = 3,
          vjust = 0,
          hjust = 0,
          check_overlap = TRUE,
          position = ggplot2::position_nudge(
            x = diff(MCB_lim) / 80,
            y = -diff(DSC_lim) / 40
          )
        )
      } else if (isTRUE(MCBDSC_repel)) {
        ggrepel::geom_text_repel(
          data = decomp_by_state$within,
          mapping = ggplot2::aes(x = .data$MCB, y = .data$DSC, label = .data$forecast, colour = .data$forecast),
          size = 3
        )
      }
    } +
    {
      if (nrow(decomp_by_state$infty)) {
        list(
          ggplot2::geom_rug(
            data = decomp_by_state$infty,
            mapping = ggplot2::aes(x = .data$MCB, y = .data$DSC, colour = .data$forecast),
            sides = "r",
            linewidth = 2
          ),
          ggplot2::geom_text(
            data = decomp_by_state$infty,
            mapping = ggplot2::aes(x = .data$x_geom_text, y = .data$DSC, label = .data$forecast, colour = .data$forecast),
            size = 3,
            hjust = 1,
            check_overlap = TRUE
          )
        )
      }
    } +
    ggplot2::scale_colour_manual(values = colour_values) +
    ggplot2::scale_x_continuous(oob = scales::oob_squish_infinite) +
    ggplot2::coord_cartesian(xlim = MCB_lim, ylim = DSC_lim) +
    ggplot2::labs(x = expression(MCB[I]), y = expression(DSC[I])) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1
      )
    )
}

# argument checks
check_plot_args <- function(decomp, n_isolines, colour_values, colour_unc, MCBDSC_repel, MCB_lim, DSC_lim) {
  # decomp
  if (!is.list(decomp) && !is.data.frame(decomp)) stop("'int' must either be a list or dataframe")
  if (is.data.frame(decomp)) {
    n <- nrow(decomp)
    if (!all(c("UNC", "DSC", "MCB", "forecast") %in% colnames(decomp)))
      stop("'decomp' must contain columns labelled 'UNC', 'DSC', 'MCB', and 'forecast'")
  } else {
    n <- length(decomp)
    if (is.null(names(decomp))) warning("no labels are given for the forecast methods")
    if (!all(sapply(decomp, function(x) all(c("UNC", "DSC", "MCB") %in% names(x)))))
      stop("the elements of 'decomp' must be named vectors with labels 'UNC', 'DSC', and 'MCB'")
  }

  # n_isolines
  if (!is.numeric(n_isolines) || length(n_isolines) > 1) stop("'n_isolines' must be a single integer")
  if ((n_isolines %% 1) != 0) stop("'n_isolines' must be an integer")

  # colour_values
  if (!is.vector(colour_values) || !is.character(colour_values)) stop("'colour_values' must be a character vector")
  if (!(length(colour_values) %in% c(1, n))) stop("'colour_values' must have length 1 or the same number of rows/elements of 'decomp'")

  # colour_unc
  if (!is.character(colour_unc) || length(colour_unc) > 1) stop("'colour_unc' must be a single character string")

  # MCBDSC_repel
  if (!is.logical(MCBDSC_repel) || length(MCBDSC_repel) > 1) stop("'MCBDSC_repel' must be a single logical")

  # MCB_lim
  if (!is.na(MCB_lim)) {
    if (!is.numeric(MCB_lim) || length(MCB_lim) > 2) stop("'MCB_lim' must be a numeric vector of length 1 or 2")
    if (length(MCB_lim) == 2) {
      if (MCB_lim[2] <= MCB_lim[1]) stop("The second value in 'MCB_lim' must be larger than the first.")
      if (MCB_lim[2] <= 0) stop(paste("The upper bound of the x-axis must be positive. Currently", MCB_lim[2]))
    } else {
      if (MCB_lim <= 0) stop(paste("The upper bound of the x-axis must be positive. Currently", MCB_lim))
    }
  }

  # DSC_lim
  if (!is.na(DSC_lim)) {
    if (!is.numeric(DSC_lim) || length(DSC_lim) > 2) stop("'DSC_lim' must be a numeric vector of length 1 or 2")
    if (length(DSC_lim) == 2) {
      if (DSC_lim[2] <= DSC_lim[1]) stop("The second value in 'DSC_lim' must be larger than the first.")
      if (DSC_lim[2] <= 0) stop(paste("The upper bound of the y-axis must be positive. Currently", DSC_lim[2]))
    } else {
      if (DSC_lim <= 0) stop(paste("The upper bound of the y-axis must be positive. Currently", DSC_lim))
    }
  }
}
