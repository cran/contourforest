#' Contour-Enhanced Binary Outcome Forest Plot
#'
#' Generates a contour-enhanced forest plot for binary outcome data (e.g., odds ratios, risk ratios),
#' with study-level effects, confidence intervals, pooled effect, prediction interval, and heterogeneity statistics.
#'
#' @param dat Data frame containing study-level binary outcome data.
#' @param measure Character. Effect measure ("OR" for odds ratio, "RR" for risk ratio, etc.).
#' @param method Character. Method for meta-analysis heterogeneity estimation (default "REML").
#' @param xlab Character. Label for the x-axis.
#' @param title Character. Plot title. If NULL, a default title is generated.
#' @param model Character. Meta-analysis model ("Random-effects" or "Fixed-effects").
#' @param estimator Character. Estimator used in the meta-analysis (default "REML").
#' @param nc_col Character. Column name for control group sample sizes.
#' @param ne_col Character. Column name for treatment group sample sizes.
#' @param event_c_col Character. Column name for number of events in control group.
#' @param event_t_col Character. Column name for number of events in treatment group.
#' @param diamond.col Color of the pooled effect polygon.
#' @param study.col Color of the study-level effect points.
#' @param CI.col Color of the study-level confidence interval lines.
#' @param Pred.Inter.col Color of the prediction interval line.
#' @param square.size Numeric. Maximum size of study-level effect squares.
#' @param contour_fill Vector of colors for contour shading levels.
#' @param text_size Numeric. Base size of plot text.
#' @param xlim Numeric vector of length 2. Limits of the x-axis.
#' @param pred Logical. Whether to show the prediction interval.
#' @param xpos List of numeric positions for text labels (EventsT, EventsC, Effect, Weight).
#' @param study_x Numeric. X-position for study names.
#' @param hetero_x Numeric. X-position for heterogeneity text.
#' @param tlim Numeric vector of length 2. Limits for truncating study confidence intervals.
#' @param truncate_PI Logical. Whether to truncate the prediction interval to `tlim`.
#' @param contour_left_min Numeric vector. Minimum x-values for left-side contour shading.
#' @param contour_left_max Numeric vector. Maximum x-values for left-side contour shading.
#' @param contour_right_min Numeric vector. Minimum x-values for right-side contour shading.
#' @param contour_right_max Numeric vector. Maximum x-values for right-side contour shading.
#'
#' @return A `ggplot2` object of the forest plot.
#'
#' @import ggplot2
#' @import metafor
#' @import dplyr
#' @import grid
#' @import stringr
#' @importFrom stats confint predict qchisq
#'
#' @export
#' @examples
#' forest_bin(
#'   dat = bcg(),
#'   measure = "OR",
#'   xlab = "Odds Ratio",
#'   title = "BCG Vaccine Meta-analysis",
#'   tlim = c(0, 2.3),
#'   contour_left_min = c(0,0.3,0.5,0.7),
#'   contour_left_max = c(0.3,0.5,0.7,1),
#'   contour_right_min = c(1,1.2,1.5,1.8),
#'   contour_right_max = c(1.2,1.5,1.8,2.5)
#' )
forest_bin <- function(
    dat,
    measure = "OR",
    method = "REML",
    xlab = "",
    title = NULL,
    model = "Random-effects",
    estimator = "REML",
    nc_col = "n_c",
    ne_col = "n_t",
    event_c_col = "events_c",
    event_t_col = "events_t",
    diamond.col = "red",
    study.col = "blue",
    CI.col = "blue",
    Pred.Inter.col = "black",
    square.size = 10,
    contour_fill = c("gray95", "gray80", "gray60", "gray40"),
    text_size = 3.5,
    xlim = c(-1.7, 3.5),
    pred = TRUE,
    xpos = list(
      EventsT = -0.9,
      EventsC = -0.3,
      Effect = 2.6,
      Weight = 3.1
    ),
    study_x = -1.8,
    hetero_x = -1.7,
    tlim = c(0,2.3),
    truncate_PI = FALSE,
    contour_left_min = c(0, 0.5, 0.67, 0.83),  # defaults
    contour_left_max = c(0.5, 0.67, 0.83, 1),
    contour_right_min = c(1, 1.2, 1.5, 2),
    contour_right_max = c(1.2, 1.5, 2, 2.5)
) {
  # --- Check required columns ---
  required_cols <- c("n_t", "n_c", "events_t", "events_c")
  missing_cols <- required_cols[!required_cols %in% names(dat)]
  if(length(missing_cols) > 0){
    stop(paste("Missing required columns in your dataset:", paste(missing_cols, collapse=", ")))
  }

  events_t <- dat[[event_t_col]]
  n_t <- dat[[ne_col]]
  events_c <- dat[[event_c_col]]
  n_c <- dat[[nc_col]]

  # Run meta-analysis
  meta_res <- metafor::rma(
    measure = measure,
    ai = events_t,
    n1i = n_t,
    ci = events_c,
    n2i = n_c,
    data = dat,
    method = method
  )

  # Prediction interval
  pi <- predict(meta_res, digits = 3, transf = exp)

  # Effect sizes, SE, CI, weights
  dat <- dat %>%
    mutate(
      yi = if (measure == "OR") {
        log((events_t / (n_t - events_t)) / (events_c / (n_c - events_c)))
      } else {
        log((events_t / n_t) / (events_c / n_c))
      },
      sei = sqrt(1/events_t + 1/(n_t - events_t) + 1/events_c + 1/(n_c - events_c)),
      weight = 100 * (1/sei^2) / sum(1/sei^2),
      ci.lb = yi - 1.96 * sei,
      ci.ub = yi + 1.96 * sei
    ) %>%
    arrange(desc(yi)) %>%
    mutate(
      y = rev(seq_len(nrow(.))),
      Effect = paste0(round(exp(yi), 2), " [", round(exp(ci.lb), 2), "-", round(exp(ci.ub), 2), "]"),
      WeightText = paste0(round(weight, 1), "%"),
      EventsT = paste0(events_t, "/", n_t),
      EventsC = paste0(events_c, "/", n_c)
    )

  # Apply truncation for selected studies
  if(!is.null(tlim)){
    log_tlim <- log(tlim + 1e-8)
    dat <- dat %>%
      mutate(
        ci.lb.trunc = pmax(ci.lb, log_tlim[1]),
        ci.ub.trunc = pmin(ci.ub, log_tlim[2]),
        arrow.left  = ci.lb < log_tlim[1],
        arrow.right = ci.ub > log_tlim[2]
      )
  } else {
    dat <- dat %>%
      mutate(
        ci.lb.trunc = ci.lb,
        ci.ub.trunc = ci.ub,
        arrow.left = FALSE,
        arrow.right = FALSE
      )
  }

  # Pooled effect polygon
  pooled <- data.frame(
    x = c(exp(meta_res$ci.lb), exp(meta_res$beta), exp(meta_res$ci.ub), exp(meta_res$beta), exp(meta_res$ci.lb)),
    y = c(0, 0.5, 0, -0.5, 0)
  )
  pooled_text <- paste0(round(exp(meta_res$beta), 2), " [", round(exp(meta_res$ci.lb), 2), "-", round(exp(meta_res$ci.ub), 2), "]")
  pi_text <- paste0("[", round(pi$pi.lb, 2), "-", round(pi$pi.ub, 2), "]")  # Prediction interval text below pooled

  # Heterogeneity
  tau_ci <- confint(meta_res)$random
  tau2 <- meta_res$tau2
  I2 <- meta_res$I2
  Q <- meta_res$QE
  df <- meta_res$k - 1
  alpha <- 0.05
  I2_lower <- max(0, (Q - qchisq(1 - alpha/2, df))/Q * 100)
  I2_upper <- min(100, (Q - qchisq(alpha/2, df))/Q * 100)

  # Determine Q-test p-value text
  Qp_text <- if(meta_res$QEp < 0.001) "<0.001" else signif(meta_res$QEp, 3)

  # Build heterogeneity text
  hetero_text <- paste0(
    model, " (", estimator, ")\n",
    "Tau^2 = ", round(tau2, 3), " [", round(tau_ci[1], 3), "-", round(tau_ci[2], 3), "]\n",
    "I^2 = ", round(I2, 1), "% [", round(I2_lower, 1), "-", round(I2_upper, 1), "%]\n",
    "Q-test p = ", Qp_text
  )
  # contour
  # Contours
  ymax_plot <- max(dat$y) + 2
  ymin_plot <- -2

  # Default values if arguments are not provided
  if(missing(contour_left_min)) contour_left_min <- c(0, 1/2, 1/1.5, 1/1.2)
  if(missing(contour_left_max)) contour_left_max <- c(1/2, 1/1.5, 1/1.2, 1)
  if(missing(contour_right_min)) contour_right_min <- c(1, 1.2, 1.5, 2)
  if(missing(contour_right_max)) contour_right_max <- c(1.2, 1.5, 2, 2.5)

  contour_left <- data.frame(
    xmin = contour_left_min,
    xmax = contour_left_max,
    ymin = ymin_plot,
    ymax = ymax_plot,
    fill = factor(c("Very Large", "Large", "Medium", "Small"),
                  levels = c("Small","Medium","Large","Very Large"))
  )

  contour_right <- data.frame(
    xmin = contour_right_min,
    xmax = contour_right_max,
    ymin = ymin_plot,
    ymax = ymax_plot,
    fill = factor(c("Small","Medium","Large","Very Large"),
                  levels = c("Small","Medium","Large","Very Large"))
  )

  contour_all <- rbind(contour_left, contour_right)



  if(is.null(title)) title <- paste0("Contour Enhanced Forest Plot (", measure, ")")

  # Plot
  p <- ggplot() +
    geom_rect(data = contour_all, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill), alpha=0.8) +
    scale_fill_manual(values=contour_fill, name="", guide=guide_legend(override.aes=list(alpha=1))) +
    geom_vline(xintercept=1, linetype="solid", alpha=0.7, color="black", linewidth=0.5) +
    geom_point(data=dat, aes(y=y, x=exp(yi), size=weight), shape=15, color=study.col, show.legend=FALSE) +
    scale_size_continuous(range=c(2,square.size)) +
    geom_errorbarh(data=dat, aes(y=y, xmin=exp(ci.lb.trunc), xmax=exp(ci.ub.trunc)), height=0, color=CI.col) +
    geom_polygon(data=pooled, aes(x=x, y=y), fill=diamond.col, alpha=0.9) +
    geom_text(data=dat, aes(x=xpos$EventsT, y=y, label=EventsT), hjust=0, size=text_size) +
    geom_text(data=dat, aes(x=xpos$EventsC, y=y, label=EventsC), hjust=0, size=text_size) +
    geom_text(data=dat, aes(x=xpos$Effect, y=y, label=Effect), hjust=0, size=text_size) +
    geom_text(data=dat, aes(x=xpos$Weight, y=y, label=WeightText), hjust=0, size=text_size) +
    annotate("text", x=study_x, y=dat$y, label=dat$Study, hjust=0, size=text_size) +
    scale_y_continuous(breaks=dat$y, labels=rep("", length(dat$y))) +
    labs(x=xlab, y=NULL, title=title) +
    theme_minimal(base_size=14) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text.y=element_blank(), legend.position="top") +
    annotate("text", x=c(xpos$EventsT,xpos$EventsC,xpos$Effect,xpos$Weight), y=rep(max(dat$y)+1,4),
             label=c("Events (T)","Events (C)",paste0(measure," (95% CI)"),"Weight (%)"), fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=study_x, y=max(dat$y)+1, label="Study", fontface="bold", hjust=0, size=text_size+0.3) +
    annotate("text", x=study_x+0.8, y=0.4, label="Pooled effect", fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=study_x+0.8, y=-0.5, label="Pred. Interval", fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=xpos$Effect, y=0.4, label=pooled_text, fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=xpos$Effect, y=-0.5, label=pi_text, fontface="bold", hjust=0, size=text_size, color="black") + # Prediction interval below pooled
    annotate("text", x=hetero_x, y=-1.1, label=hetero_text, hjust=0, fontface="italic", size=text_size)

  # Add arrows for truncated studies
  if(any(dat$arrow.left)){
    p <- p + geom_segment(data=dat %>% filter(arrow.left),
                          aes(x=exp(ci.lb.trunc), xend=exp(ci.lb.trunc)*0.95, y=y, yend=y),
                          arrow=arrow(length=unit(0.15,"cm")), color=CI.col)
  }
  if(any(dat$arrow.right)){
    p <- p + geom_segment(data=dat %>% filter(arrow.right),
                          aes(x=exp(ci.ub.trunc), xend=exp(ci.ub.trunc)*1.05, y=y, yend=y),
                          arrow=arrow(length=unit(0.15,"cm")), color=CI.col)
  }

  # Prediction interval truncation
  if(pred){
    pi.lb <- ifelse(truncate_PI & !is.null(tlim), pmax(pi$pi.lb, tlim[1]), pi$pi.lb)
    pi.ub <- ifelse(truncate_PI & !is.null(tlim), pmin(pi$pi.ub, tlim[2]), pi$pi.ub)
    p <- p + geom_errorbarh(aes(xmin=pi.lb, xmax=pi.ub, y=-0.6), height=0, size=1.9, color=Pred.Inter.col)
  }

  p <- p + coord_cartesian(xlim=xlim, ylim=c(-1.5, max(dat$y)+1.5))

  return(p)
}


