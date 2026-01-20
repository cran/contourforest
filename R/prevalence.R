#' Prevalence Forest Plot
#'
#' Generates a forest plot for prevalence data with study-level estimates, confidence intervals,
#' weights, and pooled prevalence. Supports contour shading for prevalence ranges.
#'
#' @param dat A data frame containing study-level prevalence data.
#' @param events_col Name of the column containing the number of events.
#' @param n_col Name of the column containing the total sample size.
#' @param study_col Name of the column with study labels.
#' @param xlab Label for the x-axis. Default is `"Prevalence (%)"`.
#' @param title Plot title.
#' @param xlim Numeric vector of length 2 specifying x-axis limits.
#' @param diamond.col Color of the pooled prevalence diamond.
#' @param study.col Color of the study points.
#' @param CI.col Color of the confidence intervals.
#' @param square.size Size of the study-level squares.
#' @param text_size Base text size for labels and annotations.
#' @param xpos_study X-axis position for study labels.
#' @param xpos_events X-axis position for event counts.
#' @param xpos_prev X-axis position for prevalence values.
#' @param xpos_weight X-axis position for study weights.
#' @param xpos_pooled_label X-axis position for pooled label (optional; defaults relative to pooled prevalence).
#' @param xpos_pooled_value X-axis position for pooled value (optional; defaults relative to pooled prevalence).
#' @param contour_breaks Numeric vector defining contour thresholds.
#' @param contour_colors Colors for contour shading.
#' @param legend_title Title for the contour legend.
#'
#' @return A ggplot object representing the prevalence forest plot.
#' @export
#'
#' @examples
#' dat <- data.frame(
#'   Study = c("Bazargani et al. 2015", "Nandhra et al. 2015", "Rai 2015", "Romano et al. 2011"),
#'   events = c(25, 109, 24, 9),
#'   n = c(454, 757, 380, 100)
#' )
#'
#' forest_prev(
#'   dat = dat,
#'   title = "Prevalence Forest Plot",
#'   xlim = c(-25, 80),
#'   xpos_study = -12,
#'   xpos_events = -5,
#'   xpos_prev = 40,
#'   xpos_weight = 52,
#'   contour_breaks = c(0,5,10,25,35),
#'   contour_colors = c("gray90","gray70","gray50","gray30"),
#'   legend_title = "Prevalence (%)"
#' )
forest_prev <- function(dat,
                              events_col = "events",
                              n_col = "n",
                              study_col = "Study",
                              xlab = "Prevalence (%)",
                              title = "Prevalence Forest Plot",
                              xlim = c(-25, 140),
                              diamond.col = "red",
                              study.col = "blue",
                              CI.col = "blue",
                              square.size = 7,
                              text_size = 3.5,
                              # Positions
                              xpos_study = -12,
                              xpos_events = -5,
                              xpos_prev = 50,
                              xpos_weight = 55,
                              xpos_pooled_label = NULL,
                              xpos_pooled_value = NULL,
                              # Contour thresholds
                              contour_breaks = c(0,5,10,25,50),
                              contour_colors = c("gray90","gray70","gray50","gray30"),
                              legend_title = "Prevalence"
) {

  events <- dat[[events_col]]
  n <- dat[[n_col]]
  studies <- dat[[study_col]]

  # Logit-transformed prevalence
  es <- escalc(measure = "PLO", xi = events, ni = n, data = dat)

  # Random-effects meta-analysis
  res <- rma(yi, vi, data = es, method = "REML")

  # Study-level prevalence, CI, weight
  df <- data.frame(
    Study = studies,
    Events = paste0(events, "/", n),
    Prevalence = transf.ilogit(es$yi)*100,
    CI_lb = transf.ilogit(es$yi - 1.96*sqrt(es$vi))*100,
    CI_ub = transf.ilogit(es$yi + 1.96*sqrt(es$vi))*100,
    vi = es$vi
  ) %>%
    arrange(desc(Prevalence)) %>%
    mutate(y = rev(seq_len(nrow(.))),
           weight = round(100*(1/vi)/sum(1/vi),1))

  # Pooled prevalence diamond
  pooled_val <- transf.ilogit(res$beta)*100
  pooled_ci_lb <- transf.ilogit(res$ci.lb)*100
  pooled_ci_ub <- transf.ilogit(res$ci.ub)*100
  pooled_text <- paste0(round(pooled_val,1), "% [", round(pooled_ci_lb,1), "-", round(pooled_ci_ub,1), "]")

  pooled <- data.frame(
    x = c(pooled_ci_lb, pooled_val, pooled_ci_ub, pooled_val, pooled_ci_lb),
    y = c(0, 0.5, 0, -0.5, 0)
  )

  # Contour rectangles
  contour_df <- data.frame(
    xmin = contour_breaks[-length(contour_breaks)],
    xmax = contour_breaks[-1],
    ymin = -1.5,
    ymax = max(df$y)+2,
    fill = factor(paste0(contour_breaks[-length(contour_breaks)], "-", contour_breaks[-1], "%"),
                  levels = paste0(contour_breaks[-length(contour_breaks)], "-", contour_breaks[-1], "%"))
  )

  contour_fill <- setNames(contour_colors, levels(contour_df$fill))

  # Heterogeneity
  tau2 <- res$tau2
  I2 <- res$I2
  Q <- res$QE
  df_Q <- res$k - 1
  Qp <- if(res$QEp < 0.001) "<0.001" else signif(res$QEp, 3)
  Qp <- if(res$QEp < 0.001) "<0.001" else signif(res$QEp, 3)
  hetero_text <- paste0("Tau\u00B2=", round(tau2,3),
                        ", I\u00B2=", round(I2,1), "%, \n Q(", df_Q, ")=", round(Q,2),
                        ", p=", Qp)

  # Column header y-position
  header_y <- max(df$y)+1

  # Default pooled label/value positions relative to pooled_val if not provided
  if(is.null(xpos_pooled_label)) xpos_pooled_label <- pooled_val - 15
  if(is.null(xpos_pooled_value)) xpos_pooled_value <- pooled_val + 15

  # Plot
  ggplot() +
    # Contours
    geom_rect(
      data = contour_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      alpha = 0.6 ) +
    scale_fill_manual(values = contour_fill, guide = "none") +

    # Study points (squares)
    geom_point(data = df, aes(x = Prevalence, y = y), size = square.size, color = "black", fill = study.col, shape = 22) +

    # Horizontal CI whiskers
    geom_errorbarh(data = df, aes(y = y, xmin = CI_lb, xmax = CI_ub), height = 0, color = CI.col, size = 0.8) +

    # Vertical line at pooled prevalence
    geom_segment(aes(x = pooled_val, xend = pooled_val, y = min(df$y)-0.5,
                     yend = max(df$y)+2),
                 linetype = "dotted", color = "black", size = 0.7)+

    # Pooled prevalence diamond
    geom_polygon(data = pooled, aes(x = x, y = y), fill = diamond.col, color = diamond.col, alpha = 0.8) +

    # Pooled effect label and value
    annotate("text", x = xpos_study, y = 0, label ="Pooled Effect", fontface = "bold", hjust = 1, size = text_size) +
    annotate("text", x = xpos_prev , y = 0, label = pooled_text, fontface = "bold", hjust = 0, size = text_size) +

    # Left table: Study name + Events/Total + Prevalence [95% CI]
    geom_text(data = df, aes(x = xpos_study, y = y, label = Study), hjust = 1, size = text_size) +
    geom_text(data = df, aes(x = xpos_events, y = y, label = Events), hjust = 1, size = text_size) +
    geom_text(data = df, aes(x = xpos_prev, y = y,
                             label = paste0(round(Prevalence,1), "% [",
                                            round(CI_lb,1), "-",
                                            round(CI_ub,1), "%]")),
              hjust = 0, size = text_size, color = "black") +

    # Right table: Weight only
    geom_text(data = df, aes(x = xpos_weight, y = y, label = paste0(weight, "%")), hjust = 0, size = text_size) +

    # Column headers
    annotate("text", x = xpos_study, y = header_y, label = "Study", fontface="bold", size=text_size) +
    annotate("text", x = xpos_events, y = header_y, label = "Events/Total", fontface="bold", size=text_size) +
    annotate("text", x = xpos_prev, y = header_y, label = "Prevalence \n [95% CI]", fontface="bold", size=text_size) +
    annotate("text", x = xpos_weight, y = header_y, label = "Weight", fontface="bold", size=text_size) +

    # y-axis
    scale_y_continuous(breaks = df$y, labels = rep("", nrow(df))) +

    # x-axis
    scale_x_continuous(limits = xlim) +

    # Labels and theme
    labs(x = xlab, y = NULL, title = title) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank()) +

    # Heterogeneity stats
    annotate("text", x = xpos_study-8, y = -1.5,
             label = hetero_text, size = text_size,
             fontface = "italic", hjust = 0)
}


