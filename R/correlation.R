#' Forest Plot for Correlation Coefficients
#'
#' Creates a contour-enhanced forest plot for correlation coefficients from multiple studies.
#' Applies Fisher's Z transformation for meta-analysis and displays study-level correlations
#' with 95% confidence intervals, pooled correlation, weights, and heterogeneity statistics.
#'
#' @param dat Data frame containing study-level correlation data.
#' @param r_col Character. Name of the column containing correlation coefficients. Default is `"r"`.
#' @param n_col Character. Name of the column containing sample sizes. Default is `"n"`.
#' @param study_col Character. Name of the column containing study labels. Default is `"Study"`.
#' @param xlab Character. Label for the x-axis. Default is `"Correlation (r)"`.
#' @param title Character. Plot title. Default is `"Correlation Forest Plot"`.
#' @param xlim Numeric vector of length 2. Limits for the x-axis. Default is `c(-2.5,1.5)`.
#' @param diamond.col Color for the pooled effect diamond. Default is `"red"`.
#' @param study.col Color for the study points. Default is `"blue"`.
#' @param CI.col Color for the horizontal confidence intervals. Default is `"blue"`.
#' @param square.size Numeric. Size of the squares representing study effect sizes. Default is `6`.
#' @param text_size Numeric. Size of text annotations. Default is `3.5`.
#' @param xpos_study Numeric. X-position for study labels. Default is `-1.3`.
#' @param xpos_n Numeric. X-position for sample size labels. Default is `-1.1`.
#' @param xpos_ci Numeric. X-position for correlation (95% CI) labels. Default is `1.2`.
#' @param xpos_weight Numeric. X-position for study weight labels. Default is `1.45`.
#' @param xpos_pooled_label Numeric. X-position for the pooled effect label. Default is calculated automatically.
#' @param xpos_pooled_value Numeric. X-position for the pooled effect value. Default is calculated automatically.
#' @param contour_breaks Numeric vector. Breakpoints for contour shading along the x-axis. Default is `c(-1, -0.5, -0.3,-0.1, 0, 0.1,  0.3, 0.5, 1)`.
#' @param contour_colors Character vector. Colors corresponding to contour breaks. Default is `c("gray70","gray50","gray30", "gray10", "gray70", "gray50","gray30","gray10")`.
#'
#' @importFrom stats setNames
#' @return A `ggplot` object representing the correlation forest plot.
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- data.frame(
#'   Study = c("Smith 2010","Jones 2012","Lee 2015","Kim 2018"),
#'   r = c(0.2, 0.35, -0.1, 0.5),
#'   n = c(50, 120, 80, 60)
#' )
#' forest_corr(dat)
#' }
forest_corr <- function(dat,
                              r_col = "r",
                              n_col = "n",
                              study_col = "Study",
                              xlab = "Correlation (r)",
                              title = "Correlation Forest Plot",
                              xlim = c(-2.5,1.5),
                              diamond.col = "red",
                              study.col = "blue",
                              CI.col = "blue",
                              square.size = 6,
                              text_size = 3.5,
                              # Positions
                              xpos_study = -1.3,
                              xpos_n = -1.1,
                              xpos_ci = 1.2,
                              xpos_weight = 1.45,
                              xpos_pooled_label = NULL,
                              xpos_pooled_value = NULL,
                              # Contour thresholds
                              contour_breaks = c(-1, -0.5, -0.3,-0.1, 0, 0.1,  0.3, 0.5, 1),
              contour_colors = c("gray70","gray50","gray30", "gray10", "gray70", "gray50","gray30","gray10")
) {

  r <- dat[[r_col]]
  n <- dat[[n_col]]
  studies <- dat[[study_col]]

  # Fisher's Z transform
  z <- atanh(r)
  se <- 1 / sqrt(n - 3)

  # Random-effects meta-analysis
  res <- rma(yi = z, sei = se, method = "REML")

  # Study-level correlation and CI
  df <- data.frame(
    Study = studies,
    r = r,
    n = n,
    z = z,
    se = se
  ) %>%
    mutate(
      ci_lb_z = z - 1.96*se,
      ci_ub_z = z + 1.96*se,
      ci_lb = tanh(ci_lb_z),
      ci_ub = tanh(ci_ub_z)
    ) %>%
    arrange(desc(r)) %>%
    mutate(
      y = rev(seq_len(nrow(.))),
      weight = round(100*(1/se^2)/sum(1/se^2),1)
    )

  # Pooled correlation
  pooled_val <- tanh(res$beta)
  pooled_ci_lb <- tanh(res$ci.lb)
  pooled_ci_ub <- tanh(res$ci.ub)
  pooled_text <- paste0(round(pooled_val,2), " [", round(pooled_ci_lb,2), "-", round(pooled_ci_ub,2), "]")

  pooled <- data.frame(
    x = c(pooled_ci_lb, pooled_val, pooled_ci_ub, pooled_val, pooled_ci_lb),
    y = c(0,0.5,0,-0.5,0)
  )

  # Contour rectangles
  contour_df <- data.frame(
    xmin = contour_breaks[-length(contour_breaks)],
    xmax = contour_breaks[-1],
    ymin = -1,
    ymax = max(df$y)+1,
    fill = factor(paste0(contour_breaks[-length(contour_breaks)], "-", contour_breaks[-1]))
  )

  contour_fill <- setNames(contour_colors, levels(contour_df$fill))

  # Default pooled label/value positions
  if(is.null(xpos_pooled_label)) xpos_pooled_label <- pooled_val - 0.2
  if(is.null(xpos_pooled_value)) xpos_pooled_value <- pooled_val + 0.2

  # Heterogeneity
  tau2 <- res$tau2
  I2 <- res$I2
  Q <- res$QE
  df_Q <- res$k - 1
  Qp <- if(res$QEp < 0.001) "<0.001" else signif(res$QEp, 3)
  hetero_text <- paste0("Tau\u00B2=", round(tau2,3),
                        ", I\u00B2=", round(I2,1), "%, \n Q(", df_Q, ")=", round(Q,2))

  # Column header y-position
  header_y <- max(df$y)+0.7

  # Plot
  ggplot() +
    # Contours
    geom_rect(data = contour_df, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill), alpha=0.5) +
    scale_fill_manual(values = contour_fill, guide = "none") +

    # Study points
    geom_point(data = df, aes(x=r, y=y), shape=22, size=square.size, color="black", fill=study.col) +

    # Horizontal CI
    geom_errorbarh(data=df, aes(y=y, xmin=ci_lb, xmax=ci_ub), height=0, color=CI.col, size=0.8) +

    # Pooled diamond
    geom_polygon(data=pooled, aes(x=x, y=y), fill=diamond.col, color=diamond.col, alpha=0.8) +

    # Vertical dotted line at pooled
    geom_segment(aes(x=pooled_val, xend=pooled_val, y=min(df$y)-0.5, yend=max(df$y)+0.5),
                 linetype="dotted", color="black", size=0.7) +

    # Left table: Study + n
    geom_text(data=df, aes(x=xpos_study, y=y, label=Study), hjust=1, size=text_size) +
    geom_text(data=df, aes(x=xpos_n, y=y, label=paste0(n)), hjust=1, size=text_size) +

    # Right table: r (95% CI) + Weight
    geom_text(data=df, aes(x=xpos_ci, y=y, label=paste0(round(r,2), " [", round(ci_lb,2), "-", round(ci_ub,2), "]")), hjust=0, size=text_size) +
    geom_text(data=df, aes(x=xpos_weight, y=y, label=paste0(weight,"%")), hjust=0, size=text_size) +

    # Column headers
    annotate("text", x=xpos_study, y=header_y, label="Study", fontface="bold", size=text_size) +
    annotate("text", x=xpos_n, y=header_y, label="n", fontface="bold", size=text_size) +
    annotate("text", x=xpos_ci+.05, y=header_y, label="Corr(95% CI)", fontface="bold", size=text_size) +
    annotate("text", x=xpos_weight, y=header_y, label="Weight", fontface="bold", size=text_size) +

    # Pooled label/value
    annotate("text", x=xpos_study, y=0, label="Pooled Effect", fontface="bold", hjust=1, size=text_size) +
    annotate("text", x=xpos_ci, y=0, label=pooled_text, fontface="bold", hjust=0, size=text_size) +

    # Heterogeneity
    annotate("text", x=xpos_study-1, y=min(df$y)-2, label=hetero_text, fontface="italic", size=text_size, hjust=0) +

    scale_y_continuous(breaks=df$y, labels=rep("", nrow(df))) +
    scale_x_continuous(limits=xlim) +
    labs(x=xlab, y=NULL, title=title) +
    theme_minimal(base_size=14) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y=element_blank())
}


