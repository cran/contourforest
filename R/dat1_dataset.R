#' Hypothetical dataset: GAD-7 Anxiety Scores by Socioeconomic Status
#'
#' GAD-7 scores (mean ± SD) for high SES (control) vs low SES (treatment) groups across 6 studies.
#'
#' @format A data frame with 6 rows and 8 variables:
#' \describe{
#'   \item{Study}{Study name or year}
#'   \item{mean_c}{Mean GAD-7 score in high SES (control) group}
#'   \item{sd_c}{Standard deviation of GAD-7 scores in high SES (control) group}
#'   \item{n_c}{Sample size in high SES (control) group}
#'   \item{mean_t}{Mean GAD-7 score in low SES (treatment) group}
#'   \item{sd_t}{Standard deviation of GAD-7 scores in low SES (treatment) group}
#'   \item{n_t}{Sample size in low SES (treatment) group}
#'   \item{subgroup}{Geographical subgroup of the study (e.g., Asia, Europe)}
#' }
#' @source Hypothetical data
#' @export
dat1 <- data.frame(
  Study   = c("Shell 2011", "Gill 2015", "Smith 2010", "Khan 2024", "Eliade 2021", "Manfredi 2009"),
  mean_c  = c(5.2, 4.8, 6.1, 5.5, 4.9, 5.8),  # High SES (control)
  sd_c    = c(1.8, 2.0, 2.1, 1.9, 2.2, 1.7),
  n_c     = c(50, 60, 45, 55, 50, 40),
  mean_t  = c(7.4, 5.9, 8.0, 3.8, 4.2, 8.1),  # Low SES (treatment)
  sd_t    = c(2.0, 1.9, 2.2, 2.1, 2.0, 2.3),
  n_t     = c(50, 60, 45, 55, 50, 40),
  subgroup= c("Asia", "Europe","Asia", "Europe","Asia", "Europe")
)
