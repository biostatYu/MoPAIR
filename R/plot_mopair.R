plot_mopair <- function(x, type = c("modules", "calibration", "leading_edge"), top = 20L, ...) {
  if (!inherits(x, "mopair_fit")) stop("x must be a mopair_fit object.", call. = FALSE)
  type <- match.arg(type)
  if (length(top) != 1L || !is.finite(top) || top < 1L) stop("top must be positive.", call. = FALSE)
  top <- as.integer(top)
  if (type == "modules") {
    d <- x$module_results[is.finite(x$module_results$module_p), , drop = FALSE]
    d <- utils::head(d[order(d$module_p), , drop = FALSE], top)
    if (!nrow(d)) { graphics::plot.new(); graphics::text(0.5, 0.5, "No available module evidence"); return(invisible(x)) }
    values <- -log10(pmax(d$module_p, .Machine$double.xmin))
    graphics::barplot(rev(values), names.arg = rev(d$module), horiz = TRUE, las = 1,
      col = ifelse(rev(d$selected_bh), "#B24C4C", "#6F8798"), border = NA,
      xlab = expression(-log[10](module~p)), ...)
  } else if (type == "calibration") {
    if (is.null(x$calibration)) stop("This fit has no sensitivity calibration.", call. = FALSE)
    d <- x$calibration$curve
    graphics::plot(d$observed_discoveries, d$empirical_fdr, type = "b", pch = 16,
      col = "#476F86", xlab = "Observed discoveries", ylab = "Estimated FDR", ...)
    graphics::abline(h = x$metadata$alpha, lty = 2, col = "#B24C4C")
  } else {
    d <- x$mediator_results[is.finite(x$mediator_results$composite_p), , drop = FALSE]
    d <- utils::head(d[order(d$composite_p), , drop = FALSE], top)
    if (!nrow(d)) { graphics::plot.new(); graphics::text(0.5, 0.5, "No available member evidence"); return(invisible(x)) }
    values <- -log10(pmax(d$composite_p, .Machine$double.xmin))
    graphics::dotchart(rev(values), labels = rev(d$mediator), color = ifelse(rev(d$leading_edge), "#B24C4C", "#6F8798"),
      xlab = expression(-log[10](composite~statistic)), ...)
  }
  invisible(x)
}

plot.mopair_fit <- function(x, ...) plot_mopair(x, ...)
