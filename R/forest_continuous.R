#' Contour-enhanced Forest Plot for Continuous Outcomes
#'
#' Creates a forest plot for continuous outcomes (MD or SMD) with optional contour shading,
#' prediction interval, and annotated study-level means, effect sizes, and weights.
#'
#' @param dat Data frame containing study-level data. Must include treatment and control means, SDs, sample sizes, and study labels.
#' @param measure Effect size measure: `"MD"` (mean difference) or `"SMD"` (standardized mean difference). Default is `"SMD"`.
#' @param method Meta-analysis method for `rma()` (e.g., `"REML"`). Default is `"REML"`.
#' @param sort  Logical. If `TRUE`, the studies will be sorted by effect size before plotting.
#' @param xlab Label for the x-axis.
#' @param title Plot title. If `NULL`, a default title including measure is used.
#' @param model Model description for heterogeneity annotation. Default `"Random-effects"`.
#' @param estimator Estimator for heterogeneity. Default `"REML"`.
#' @param m_c_col Column name for control group means. Default `"mean_c"`.
#' @param sd_c_col Column name for control group SDs. Default `"sd_c"`.
#' @param n_c_col Column name for control group sample sizes. Default `"n_c"`.
#' @param m_t_col Column name for treatment group means. Default `"mean_t"`.
#' @param sd_t_col Column name for treatment group SDs. Default `"sd_t"`.
#' @param n_t_col Column name for treatment group sample sizes. Default `"n_t"`.
#' @param diamond.col Color for the pooled effect diamond. Default `"red"`.
#' @param study.col Color for study effect points. Default `"blue"`.
#' @param CI.col Color for study confidence intervals. Default `"blue"`.
#' @param Pred.Inter.col Color for prediction interval. Default `"black"`.
#' @param square.size Size of study points. Default `10`.
#' @param contour_fill Vector of four colors for contour shading. Default `c("gray95","gray80","gray60","gray40")`.
#' @param text_size Size of annotated text. Default `3.5`.
#' @param contour_left_min Numeric vector specifying the left-side minimum x-axis
#'   boundaries for contour shading bands (values less than the null effect).
#' @param contour_left_max Numeric vector specifying the left-side maximum x-axis
#'   boundaries for contour shading bands (values less than or equal to the null effect).
#' @param contour_right_min Numeric vector specifying the right-side minimum x-axis
#'   boundaries for contour shading bands (values greater than or equal to the null effect).
#' @param contour_right_max Numeric vector specifying the right-side maximum x-axis
#'   boundaries for contour shading bands (values greater than the null effect).
#' @param pred Logical; whether to show prediction interval. Default `TRUE`.
#' @param study_x X-position for study labels. Default computed automatically.
#' @param treatment_x X-position for treatment means. Default computed automatically.
#' @param control_x X-position for control means. Default computed automatically.
#' @param effect_x X-position for effect sizes. Default computed automatically.
#' @param weight_x X-position for weights. Default computed automatically.
#' @param PredInt_x X-position for prediction interval label. Default computed automatically.
#' @param xlim X-axis limits. Default computed automatically.
#' @param hetero_x X-position for heterogeneity annotation. Default `-8`.
#'
#' @return A `ggplot2` object representing the forest plot.
#' @export
#
#' @examples
#' forest_cont(
#'   dat1,
#'   measure = "MD",
#'   xlab = "Mean Difference",
#'   study_x = -9,
#'   sort = "effect",
#'   hetero_x = -12,
#'   treatment_x = -7,
#'   control_x = -5,
#'   effect_x = 5.5,
#'   weight_x = 10,
#'   PredInt_x = 7
#' )
#'
#' forest_cont(
#'   dat1,
#'   measure = "SMD",
#'   xlab = "Standardized Mean Difference",
#'   hetero_x = -9.9,
#'   study_x = -7,
#'   sort = "effect",
#'   treatment_x = -5,
#'   control_x = -3,
#'   effect_x = 2.5,
#'   weight_x = 4,
#'   PredInt_x = 4
#' )
#' forest_cont(
#' dat1,
#' study.col = "darkgreen",
#' CI.col = "black",
#' diamond.col = "red",
#' Pred.Inter.col = "black",
#' measure = "SMD",
#' sort = "effect",
#' xlab = "Standardized Mean Difference",
#' contour_fill = c("gray90","gray70","gray50", "gray30"),
#' hetero_x = -9.9,
#' study_x = -7,
#' square.size = 9,
#' treatment_x = -5,
#' control_x = -3.2,
#' text_size = 4,
#' effect_x = 2.5,
#' weight_x = 5.8,
#' PredInt_x = 3
#' )
forest_cont <- function(dat,
                        measure = "SMD",
                        method = "REML",
                        sort = c("effect", "none"),
                        xlab = "",
                        title = NULL,
                        model = "Random-effects",
                        estimator = "REML",
                        m_c_col = "mean_c",
                        sd_c_col = "sd_c",
                        n_c_col = "n_c",
                        m_t_col = "mean_t",
                        sd_t_col = "sd_t",
                        n_t_col = "n_t",
                        diamond.col = "red",
                        study.col = "blue",
                        CI.col = "blue",
                        Pred.Inter.col = "black",
                        # Contour limits
                        contour_left_min  = c(-0.8, -0.5, -0.2, 0),
                        contour_left_max  = c(-0.5, -0.2,  0,   0.2),
                        contour_right_min = c(0, 0.2, 0.5, 0.8),
                        contour_right_max = c(0.2, 0.5, 0.8, 1),
                        square.size = 10,
                        contour_fill = c("gray95","gray80","gray60","gray40"),
                        text_size = 3.5,
                        pred = TRUE,
                        study_x = NULL,
                        treatment_x = NULL,
                        control_x = NULL,
                        effect_x = NULL,
                        weight_x = NULL,
                        PredInt_x = NULL,
                        xlim = NULL,
                        hetero_x = -8) {

  # -------------------------
  # Compute effect size
  # -------------------------
  escalc_res <- metafor::escalc(measure=measure,
                                m1i=dat[[m_t_col]], sd1i=dat[[sd_t_col]], n1i=dat[[n_t_col]],
                                m2i=dat[[m_c_col]], sd2i=dat[[sd_c_col]], n2i=dat[[n_c_col]],
                                data=dat)

  meta_res <- metafor::rma(yi, vi, data=escalc_res, method=method)
  pi <- predict(meta_res, digits=3)

  # -------------------------
  # Prepare data for plotting
  # -------------------------
  dat_plot <- escalc_res %>%
    mutate(
      weight = 100 * (1/vi)/sum(1/vi),
      ci.lb = yi - 1.96*sqrt(vi),
      ci.ub = yi + 1.96*sqrt(vi),
      Effect = sprintf("%.2f [%.2f - %.2f]", yi, ci.lb, ci.ub),
      WeightText = sprintf("%.1f%%", weight),
      MeanT = sprintf("%.2f [%.2f]", dat[[m_t_col]], dat[[sd_t_col]]),
      MeanC = sprintf("%.2f [%.2f]", dat[[m_c_col]], dat[[sd_c_col]])
    ) %>%
    {
      if (sort == "effect") arrange(., desc(yi)) else .
    } %>%
    mutate(y = rev(seq_len(nrow(.))))


  pooled <- data.frame(
    x = c(meta_res$ci.lb, meta_res$beta, meta_res$ci.ub, meta_res$beta, meta_res$ci.lb),
    y = c(0, 0.5, 0, -0.5, 0)
  )
  pooled_text <- sprintf("%.2f [%.2f - %.2f]", meta_res$beta, meta_res$ci.lb, meta_res$ci.ub)

  tau_ci <- confint(meta_res)$random
  tau2 <- meta_res$tau2
  I2 <- meta_res$I2
  Q <- meta_res$QE
  df <- meta_res$k - 1
  I2_lb <- max(0, (Q/df - 1)/Q * 100)
  I2_ub <- min(100, ((Q * qchisq(0.975, df)/df) - 1)/Q * 100)
  hetero_text <- paste0(model, " (", estimator, ")\n",
                        "Tau^2 = ", round(tau2,3), " [", round(tau_ci[1],3), " - ", round(tau_ci[2],3), "]\n",
                        "I^2 = ", round(I2,1), "% [", round(I2_lb,1), " - ", round(I2_ub,1), "%]\n",
                        "Q-test p = ", signif(meta_res$QEp,3))

  # -------------------------
  # Default column positions
  # -------------------------
  effect_min <- pi$pi.lb
  effect_max <- pi$pi.ub

  if(is.null(study_x)) study_x <- effect_min - 4
  if(is.null(treatment_x)) treatment_x <- effect_min - 3
  if(is.null(control_x)) control_x <- effect_min - 1.5
  if(is.null(effect_x)) effect_x <- effect_max + 0.3
  if(is.null(weight_x)) weight_x <- effect_max + 3
  if(is.null(PredInt_x)) PredInt_x <- effect_max + 0.9

  if(is.null(hetero_x)) hetero_x <- min(study_x, effect_min) - 0.5
  if(is.null(xlim)) xlim <- c(min(study_x-1, hetero_x-0.5), weight_x + 2)

  # -------------------------
  # Contour shading
  # -------------------------
  # -------------------------
  # Contour shading (user-defined)
  # -------------------------
  ymax_plot <- max(dat_plot$y) + 2
  ymin_plot <- -2

  # Extend first & last bands automatically
  contour_left_min[1]  <- min(contour_left_min[1], effect_min)
  contour_right_max[length(contour_right_max)] <-
    max(contour_right_max[length(contour_right_max)], effect_max)

  contour_left <- data.frame(
    xmin = contour_left_min,
    xmax = contour_left_max,
    ymin = ymin_plot,
    ymax = ymax_plot,
    fill = factor(
      rev(c("Small","Medium","Large","Very Large")),
      levels = c("Small","Medium","Large","Very Large")
    )
  )

  contour_right <- data.frame(
    xmin = contour_right_min,
    xmax = contour_right_max,
    ymin = ymin_plot,
    ymax = ymax_plot,
    fill = factor(
      c("Small","Medium","Large","Very Large"),
      levels = c("Small","Medium","Large","Very Large")
    )
  )

  contour_all <- rbind(contour_left, contour_right)


  pred_y <- -0.6

  # -------------------------
  # Build plot
  # -------------------------
  p <- ggplot() +
    geom_rect(data = contour_all, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.7) +
    scale_fill_manual(values = contour_fill, name = "", guide = guide_legend(override.aes = list(alpha = 1))) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    geom_errorbarh(data = dat_plot, aes(y = y, xmin = ci.lb, xmax = ci.ub), height = 0, color = CI.col) +
    geom_point(data = dat_plot, aes(y = y, x = yi, size = weight), shape = 15, color = study.col, show.legend = FALSE) +
    scale_size_continuous(range = c(1, square.size))+
    geom_polygon(data = pooled, aes(x = x, y = y), fill = diamond.col, alpha = 0.9)

  if(pred){
    p <- p +
      geom_errorbar(data = data.frame(y=pred_y, xmin=pi$pi.lb, xmax=pi$pi.ub),
                    aes(y=y, xmin=xmin, xmax=xmax),
                    height = 0, size = 1.5, color = Pred.Inter.col) +
      annotate("text", x = PredInt_x, y = pred_y,
               label = paste0("[", round(pi$pi.lb,2), "-", round(pi$pi.ub,2), "]"),
               fontface="bold", hjust=0, size=text_size) +
      annotate("text", x = study_x+.9, y = pred_y,
               label = "Prediction Interval", fontface="bold", hjust=0, size=text_size)
  }

  # -------------------------
  # Left & Right tables
  # -------------------------
  p <- p +
    annotate("text", x=study_x, y=dat_plot$y, label=dat_plot$Study, hjust=0, size=text_size, family="mono") +
    annotate("text", x=treatment_x, y=dat_plot$y, label=dat_plot$MeanT, hjust=0, size=text_size, family="mono") +
    annotate("text", x=control_x, y=dat_plot$y, label=dat_plot$MeanC, hjust=0, size=text_size, family="mono") +
    annotate("text", x=effect_x, y=dat_plot$y, label=dat_plot$Effect, hjust=0, size=text_size, family="mono") +
    annotate("text", x=weight_x, y=dat_plot$y, label=dat_plot$WeightText, hjust=0, size=text_size, family="mono") +

    # Headers
    annotate("text", x=study_x-0.3, y=max(dat_plot$y)+1, label="Study", fontface="bold", hjust=0, size=text_size+.3) +
    annotate("text", x=treatment_x, y=max(dat_plot$y)+1, label="Treatment \nmean[sd]", fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=control_x, y=max(dat_plot$y)+1, label="Control \nmean[sd]", fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=effect_x, y=max(dat_plot$y)+1, label=paste0(measure," (95% CI)"), fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=weight_x, y=max(dat_plot$y)+1, label="Weight (%)", fontface="bold", hjust=0, size=text_size) +

    # Pooled effect
    annotate("text", x=study_x+.9, y=0.3, label="Pooled effect", fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=effect_x, y=0.4, label=pooled_text, fontface="bold", hjust=0, size=text_size) +

    # Heterogeneity
    annotate("text", x=hetero_x, y=-1.1, label=hetero_text, hjust=0, fontface="italic", size=text_size) +

    scale_y_continuous(breaks = dat_plot$y, labels = rep("", length(dat_plot$y))) +
    labs(x=xlab, y=NULL, title=ifelse(is.null(title), paste0("Contour Enhanced Forest Plot (", measure, ")"), title)) +
    theme_minimal(base_size=14, base_family="mono") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "top") +
    coord_cartesian(xlim=xlim, ylim=c(pred_y-0.5, max(dat_plot$y)+1.5))

  return(p)
}




