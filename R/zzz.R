#' @title Message for the User
#'
#' @description Prints out a friendly reminder message to the user.
#'
#' @usage NULL
#' @return NULL
#'
.onAttach <- function(...) {
  ver <- utils::packageVersion("bnmonitor")
  packageStartupMessage("This is bnmonitor version ", ver)
  packageStartupMessage("")
  packageStartupMessage("- If you are using bnmonitor, remember to cite:")
  packageStartupMessage("")
  packageStartupMessage("Leonelli, M., Ramanathan, R., & Wilkerson, R. L. (2023). Sensitivity and robustness analysis in Bayesian networks with the bnmonitor R package. Knowledge-Based Systems, 278, 110882.")
  packageStartupMessage("")
}
