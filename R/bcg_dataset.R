#' BCG Vaccine Trials Dataset
#'
#' This dataset contains summary data from 13 BCG vaccine studies, including the number of events
#' in treatment and control groups, total sample sizes, and subgroup classification.
#'
#' @return A data frame with 13 rows and 6 variables:
#' \describe{
#'   \item{Study}{Character. Name and year of the study.}
#'   \item{events_t}{Numeric. Number of events (e.g., TB cases) in the treatment (BCG) group.}
#'   \item{n_t}{Numeric. Total number of participants in the treatment group.}
#'   \item{events_c}{Numeric. Number of events in the control group.}
#'   \item{n_c}{Numeric. Total number of participants in the control group.}
#'   \item{subgroup}{Character. Type of study or subgroup classification ("random", "alternate", or "systematic").}
#' }
#' @export
#'
bcg <- function() {
  data.frame(
    Study    = c(
      "Aronson 1948", "Ferguson & Simes 1949", "Rosenthal et al 1960",
      "Hart & Sutherland 1977", "Frimodt-Moller et al 1973", "Stein & Aronson 1953",
      "Vandiviere et al 1973", "TPT Madras 1980", "Coetzee & Berjak 1968",
      "Rosenthal et al 1961", "Comstock et al 1974", "Comstock & Webster 1969",
      "Comstock et al 1976"
    ),
    events_t = c(4,6,3,62,33,180,8,505,29,17,186,5,27),
    n_t      = c(119,300,228,13536,5036,1361,2537,87886,7470,1699,50448,2493,16886),
    events_c = c(11,29,11,248,47,372,10,499,45,65,141,3,29),
    n_c      = c(128,274,209,12619,5761,1079,619,87892,7232,1600,27197,2338,17825),
    subgroup = c(
      "random","random","random","random","alternate","alternate","random",
      "random","random","systematic","systematic","systematic","systematic"
    ),
    stringsAsFactors = FALSE
  )
}
