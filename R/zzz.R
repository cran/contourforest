# zzz.R
# Suppress NOTES about global variables in NSE
utils::globalVariables(
  c(
    ".", "yi", "vi", "ci.lb","ci.lb.trunc" , "ci.ub.trunc", "sei",  "ci.ub", "weight", "x", "y",
    "xmin", "xmax", "ymin", "ymax", "fill",
    "EventsT", "EventsC", "Effect", "WeightText",
    "arrow.left", "arrow.right",
    "events_t", "events_c", "n_t", "n_c",
    "y_label", "label"
  )
)
