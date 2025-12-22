#' Subgroup Forest Plot for Binary Outcome Meta-analysis
#'
#' Generates a contour-enhanced forest plot for binary outcome data, with optional subgroup analysis.
#' This function calculates risk ratios (or other measures) and displays study-specific effect sizes,
#' pooled subgroup effects, overall pooled effect, and heterogeneity statistics (I^2, tau^2).
#'
#' @param dat A data frame containing study-level data. Must include columns for treatment and control events and sample sizes, and optionally a `subgroup` column.
#' @param subgroup Column name in `dat` specifying subgroup membership (character or factor). Defaults to NULL (no subgroups).
#' @param measure Effect measure to use. Default is "RR" (risk ratio).
#' @param method Method for random-effects meta-analysis. Default is "REML".
#' @param xlab Label for the x-axis. Default is "Risk Ratio (RR)".
#' @param title Plot title. Default is "Subgroup Forest Plot".
#' @param diamond.col Color for subgroup pooled effect diamonds. Default is "red".
#' @param overall.col Color for overall pooled effect diamond. Default is "darkgreen".
#' @param study.col Color for individual study points. Default is "blue".
#' @param nc_col Character. Column name for control group sample sizes.
#' @param ne_col Character. Column name for treatment group sample sizes.
#' @param event_c_col Character. Column name for number of events in control group.
#' @param event_t_col Character. Column name for number of events in treatment group.
#' @param CI.col Color for study confidence interval bars. Default is "blue".
#' @param Pred.Inter.col Color for prediction interval bars. Default is "black".
#' @param square.size Maximum size of study points. Default is 8.
#' @param Pred.Int.size Thickness of prediction interval line. Default is 2.
#' @param xlim Numeric vector of length 2 giving x-axis limits. Default is c(-2, 3.5).
#' @param tlim Numeric vector of length 2 for truncating study CIs. Default is c(0, 2.3).
#' @param text_size Base text size for labels. Default is 3.5.
#' @param xpos List of x-axis positions for EventsT, EventsC, Effect, and Weight labels. Default is list(EventsT=-1, EventsC=-0.4, Effect=2.5, Weight=3.2).
#' @param study_x X position for study names. Default is -1.8.
#' @param val_x X position for pooled effect labels. Default is 2.6.
#' @param contour_left_min Numeric vector defining left contour minimum values. Default is c(0, 0.5, 0.67, 0.83).
#' @param contour_left_max Numeric vector defining left contour maximum values. Default is c(0.5, 0.67, 0.83, 1).
#' @param contour_right_min Numeric vector defining right contour minimum values. Default is c(1, 1.2, 1.5, 2).
#' @param contour_right_max Numeric vector defining right contour maximum values. Default is c(1.2, 1.5, 2, 2.5).
#' @param pred Logical indicating whether to show the prediction interval. Default is TRUE.
#'
#' @return A `ggplot` object representing the contour-enhanced subgroup forest plot.
#'
#' @export
#'
#' @examples
#' # Load example dataset
#' data <- bcg()
#'
#' # Generate subgroup forest plot
#' forest_bin_subgroup(
#'   dat = bcg(),
#'   tlim = c(0, 2.3),
#'   contour_left_min  = c(0, 0.3, 0.5, 0.7),
#'   contour_left_max  = c(0.3, 0.5, 0.7, 1),
#'   contour_right_min = c(1, 1.2, 1.5, 1.8),
#'   contour_right_max = c(1.2, 1.5, 1.8, 2.4)
#' )
forest_bin_subgroup <- function(
    dat,
    subgroup = NULL,
    measure = "RR",
    method = "REML",
    nc_col = "n_c",
    ne_col = "n_t",
    event_c_col = "events_c",
    event_t_col = "events_t",
    xlab = "Risk Ratio (RR)",
    title = "Subgroup Forest Plot",
    diamond.col = "red",
    overall.col = "darkgreen",
    study.col = "blue",
    CI.col = "blue",
    Pred.Inter.col = "black",
    square.size = 8,
    Pred.Int.size = 2,
    xlim = c(-2, 3.5),
    tlim = c(0, 2.3),
    text_size = 3.5,
    xpos = list(EventsT = -1, EventsC = -0.4, Effect = 2.5, Weight = 3.2),
    study_x = -1.8,
    val_x = 2.6,
    contour_left_min = c(0, 0.5, 0.67, 0.83),  # defaults
    contour_left_max = c(0.5, 0.67, 0.83, 1),
    contour_right_min = c(1, 1.2, 1.5, 2),
    contour_right_max = c(1.2, 1.5, 2, 2.5),
    pred = TRUE
)
{
  # --- Check required columns ---
  required_cols <- c("n_t", "n_c", "events_t", "events_c")
  missing_cols <- required_cols[!required_cols %in% names(dat)]
  if(length(missing_cols) > 0){
    stop(paste("Missing required columns in your dataset:", paste(missing_cols, collapse=", ")))
  }

  dat <- dat %>%
    mutate(
      yi = log((events_t / n_t) / (events_c / n_c)),
      sei = sqrt(1/events_t - 1/n_t + 1/events_c - 1/n_c),
      weight = 100 * (1/sei^2) / sum(1/sei^2),
      ci.lb = yi - 1.96 * sei,       # original lower CI
      ci.ub = yi + 1.96 * sei        # original upper CI
    )

  # Create truncated CIs for plotting without changing original CIs
  if (!is.null(tlim)) {
    log_tlim <- log(tlim + 1e-8)  # avoid log(0)
    dat <- dat %>%
      mutate(
        ci.lb.trunc = pmax(ci.lb, log_tlim[1]),
        ci.ub.trunc = pmin(ci.ub, log_tlim[2]),
        arrow.left = ci.lb < log_tlim[1],
        arrow.right = ci.ub > log_tlim[2]
      )
  } else {
    dat <- dat %>%
      mutate(ci.lb.trunc = ci.lb, ci.ub.trunc = ci.ub,
             arrow.left = FALSE, arrow.right = FALSE)
  }



  sub_results <- dat %>%
    group_by(subgroup) %>%
    group_modify(~{
      m <- rma(ai = .x$events_t, n1i = .x$n_t, ci = .x$events_c, n2i = .x$n_c,
               data = .x, measure = measure, method = method)
      tau_ci <- confint(m, parm = "tau2")$random
      I2.se <- sqrt(2 / (m$k - 1)) * 100
      tibble(
        subgroup = unique(.x$subgroup),
        yi = round(m$b, 2),
        ci.lb = round(m$ci.lb, 2),
        ci.ub = round(m$ci.ub, 2),
        tau2 = round(m$tau2, 2),
        tau2.lb = round(tau_ci[1], 2),
        tau2.ub = round(tau_ci[2], 2),
        I2 = round(m$I2, 2),
        I2.lb = round(max(0, m$I2 - 1.96 * I2.se), 2),
        I2.ub = round(min(100, m$I2 + 1.96 * I2.se), 2),
        pval = round(m$QEp, 3)
      )
    })

  overall <- rma(ai = events_t, n1i = n_t, ci = events_c, n2i = n_c,
                 data = dat, measure = measure, method = method)
  overall_tau_ci <- confint(overall, parm = "tau2")$random
  I2.se_overall <- sqrt(2 / (overall$k - 1)) * 100
  overall_I2_lb <- round(max(0, overall$I2 - 1.96 * I2.se_overall), 2)
  overall_I2_ub <- round(min(100, overall$I2 + 1.96 * I2.se_overall), 2)

  if (pred) {
    pi <- predict(overall, digits = 3, transf = exp)
  }

  dat <- dat %>% arrange(subgroup)
  dat$y <- NA
  y_counter <- 0
  subgroup_gap <- 2

  for (g in unique(dat$subgroup)) {
    idx <- which(dat$subgroup == g)
    n <- length(idx)
    dat$y[idx] <- seq(from = y_counter + n, to = y_counter + 1, by = -1)
    y_counter <- y_counter + n + subgroup_gap
  }

  diamonds <- list()
  hetero <- data.frame()

  for (g in unique(dat$subgroup)) {
    subset_g <- dat %>% filter(subgroup == g)
    res <- sub_results %>% filter(subgroup == g)
    y_diamond <- min(subset_g$y) - 0.7
    diamonds[[g]] <- data.frame(
      x = c(exp(res$ci.lb), exp(res$yi), exp(res$ci.ub), exp(res$yi), exp(res$ci.lb)),
      y = c(y_diamond, y_diamond + 0.4, y_diamond, y_diamond - 0.4, y_diamond)
    )
    hetero <- rbind(
      hetero,
      data.frame(
        subgroup = g,
        y = y_diamond - 0.8,
        label = paste0(
          "I^2 = ", res$I2, "% [", res$I2.lb, "-", res$I2.ub, "], tua^2 = ",
          res$tau2, " [", res$tau2.lb, "-", res$tau2.ub, "], p = ", res$pval
        )
      )
    )
  }

  y_overall <- -1.5
  overall_diamond <- data.frame(
    x = c(exp(overall$ci.lb), exp(round(overall$b, 2)), exp(overall$ci.ub), exp(round(overall$b, 2)), exp(overall$ci.lb)),
    y = c(y_overall, y_overall + 0.4, y_overall, y_overall - 0.4, y_overall)
  )

  overall_label <- paste0(
    "I^2 = ", round(overall$I2, 2), "% [", overall_I2_lb, "-", overall_I2_ub, "], tua^2 = ",
    round(overall$tau2, 2), " [", round(overall_tau_ci[1], 2), "-", round(overall_tau_ci[2], 2), "]"
  )

  model_label <- paste0("Random effects (", method, ")")

  # contour
  ymax_plot <- max(dat$y) + 2
  ymin_plot <- -3

  if (!is.null(contour_left_min) & !is.null(contour_left_max)) {
    contour_left <- data.frame(
      xmin = contour_left_min,
      xmax = contour_left_max,
      ymin = ymin_plot,
      ymax = ymax_plot,
      fill = factor(rev(c("Small", "Medium", "Large", "Very Large")), levels = c("Small", "Medium", "Large", "Very Large"))
    )
  } else {
    contour_left <- data.frame(
      xmin = c(0, 1/2, 1/1.5, 1/1.2),
      xmax = c(1/2, 1/1.5, 1/1.2, 1),
      ymin = ymin_plot,
      ymax = ymax_plot,
      fill = factor(c("Very Large", "Large", "Medium", "Small"), levels = c("Small", "Medium", "Large", "Very Large"))
    )
  }

  if (!is.null(contour_right_min) & !is.null(contour_right_max)) {
    contour_right <- data.frame(
      xmin = contour_right_min,
      xmax = contour_right_max,
      ymin = ymin_plot,
      ymax = ymax_plot,
      fill = factor(c("Small", "Medium", "Large", "Very Large"), levels = c("Small", "Medium", "Large", "Very Large"))
    )
  } else {
    contour_right <- data.frame(
      xmin = c(1, 1.2, 1.5, 2),
      xmax = c(1.2, 1.5, 2, 2.5),
      ymin = ymin_plot,
      ymax = ymax_plot,
      fill = factor(c("Small", "Medium", "Large", "Very Large"), levels = c("Small", "Medium", "Large", "Very Large"))
    )
  }

  contour_all <- rbind(contour_left, contour_right)

  # j
  group_labels <- dat %>% group_by(subgroup) %>% summarise(y_label = max(y) + 0.8)
  subgroup_labels <- sub_results %>% left_join(group_labels, by = "subgroup") %>% mutate(label = subgroup)

  p <- ggplot(dat, aes(x = exp(yi), y = y)) +
    geom_rect(
      data = contour_all,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      alpha = 0.9,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = c(Small = "gray95", Medium = "gray85", Large = "gray70", `Very Large` = "gray50"),
      name = ""
    ) +
    geom_segment(
      aes(x = 1, xend = 1, y = ymin_plot, yend = ymax_plot),
      linetype = "solid",
      color = "black",
      size = 0.8,
      inherit.aes = FALSE
    ) +
    geom_errorbarh(data = dat,
                   aes(y = y, xmin = exp(ci.lb.trunc), xmax = exp(ci.ub.trunc)),
                   height = 0, color = CI.col) +
    geom_point(aes(size = weight), shape = 15, color = study.col) +
    scale_size_continuous(range = c(2, square.size), guide = "none") +
    annotate("text", x = study_x, y = max(dat$y) + 1.5, label = "Study", hjust = 0, size = text_size + 0.3, fontface = "bold") +
    annotate("text", x = c(xpos$EventsT, xpos$EventsC, xpos$Effect, xpos$Weight),
             y = max(dat$y) + 1.5,
             label = c("Events (T)", "Events (C)", paste0(measure, " (95% CI)"), "Weight (%)"),
             fontface = "bold", hjust = 0, size = text_size + 0.3) +
    geom_text(aes(x = xpos$EventsT, label = paste0(events_t, "/", n_t)), hjust = 0, size = text_size) +
    geom_text(aes(x = xpos$EventsC, label = paste0(events_c, "/", n_c)), hjust = 0, size = text_size) +
    geom_text(aes(x = xpos$Effect, label = paste0(round(exp(yi), 2), " [", round(exp(ci.lb), 2), "-", round(exp(ci.ub), 2), "]")), hjust = 0, size = text_size) +
    geom_text(aes(x = xpos$Weight, label = paste0(round(weight, 1), "%")), hjust = 0, size = text_size) +
    annotate("text", x = study_x, y = dat$y, label = dat$Study, hjust = 0, size = text_size) +
    lapply(diamonds, function(df) geom_polygon(data = df, aes(x = x, y = y), fill = diamond.col, alpha = 0.8, inherit.aes = FALSE)) +
    geom_text(data = subgroup_labels, aes(x = study_x, y = y_label, label = paste0("Subgroup: ", label)), hjust = 0, size = text_size, fontface = "bold", inherit.aes = FALSE) +
    lapply(names(diamonds), function(g) {
      df <- diamonds[[g]]
      res <- sub_results %>% filter(subgroup == g)
      y_d <- df$y[1]
      geom_text(aes(x = val_x, y = y_d, label = paste0(round(exp(res$yi), 2), " [", round(exp(res$ci.lb), 2), "-", round(exp(res$ci.ub), 2), "]")),
                hjust = 0, size = text_size-.2, fontface = "bold", inherit.aes = FALSE)}) +
    geom_text(data = hetero, aes(x = study_x, y = y, label = label), hjust = 0, size = text_size - 0.2, fontface = "italic",
              inherit.aes = FALSE) +
    geom_polygon(data = overall_diamond, aes(x = x, y = y), fill = overall.col, alpha = 0.9, inherit.aes = FALSE) +
    annotate("text", x = study_x, y = y_overall + 0.2, label = "Overall", hjust = 0, size = text_size + 0.2, fontface = "bold") +
    annotate("text", x = study_x, y = y_overall - 1.8, label = model_label, hjust = 0, size = text_size-0.2, fontface = "bold", color = "black") +
    annotate("text", x = val_x, y = y_overall +0.2, label = paste0(round(exp(overall$b), 2), " [", round(exp(overall$ci.lb), 2), "-",
                                                                   round(exp(overall$ci.ub), 2), "]"), hjust = 0, size = text_size , fontface = "bold") +
    annotate("text", x = study_x, y = y_overall - 0.6, label = overall_label, hjust = 0, size = text_size - 0.2, fontface = "italic") +
    labs(x = xlab, y = NULL, title = title) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(), legend.position = "top") +
    coord_cartesian(xlim = xlim, ylim = c(ymin_plot, max(dat$y) + 3))

  if (pred) {
    p <- p +
      geom_errorbar(aes(xmin = pi$pi.lb, xmax = pi$pi.ub, y = y_overall - 0.8), height = 0.001, size = Pred.Int.size, color = Pred.Inter.col, inherit.aes = FALSE) +
      annotate("text", x = study_x, y = y_overall - 1.1, label = "Prediction interval",
               fontface = "bold", hjust = 0, size = text_size-0.2) +
      annotate("text", x = val_x, y = y_overall - 1,
               label = paste0("[", round(pi$pi.lb, 2), "-", round(pi$pi.ub, 2), "]"),
               fontface = "bold", color = "black", hjust = 0, size = text_size)
  }

  return(p)
}


