#' Christchurch Health and Development Study
#'
#' Simulated data and Bayesian networks from the Christchurch Health and Development Study
#'
#' @docType data
#'
#'
#'
#' @format The dataframe \code{chds} includes 500 observations randomly simulated from the \code{bn.fit} object \code{chds_bn.fit}. It has four variables:
#' \itemize{
#'    \item \bold{Social}: family's social background with levels \code{"High"} and \code{"Low"}
#'    \item \bold{Economic}: family's economic status with levels \code{"High"} and \code{"Low"}
#'    \item \bold{Events}: number of family life events with levels \code{"High"}, \code{"Average"} and \code{"Low"}
#'    \item \bold{Admission}: hospital admission of the child with levels \code{"yes"} and \code{"no"}
#'    \item \bold{statistics}: mark out of 100 for statistics
#'}
#' \code{chds_bn} is an object of class \code{bn} including the MAP Bayesian network from Barclay et al. (2013) and \code{chds_bn.fit} is an object of class \code{bn.fit} including the probabilities from the same article.
#'
#' @name chds
#'
#' @references Fergusson, D. M., Horwood, L. J., & Shannon, F. T. (1986). Social and family factors in childhood hospital admission. Journal of Epidemiology & Community Health, 40(1), 50-58.
#' @references Barclay, L. M., Hutton, J. L., & Smith, J. Q. (2013). Refining a Bayesian network using a chain event graph. International Journal of Approximate Reasoning, 54(9), 1300-1309.
#'
NULL

#' @rdname chds
"chds"

#' @rdname chds
"chds_bn"

#' @rdname chds
"chds_bn.fit"


