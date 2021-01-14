#' Bayesian network on travel survey
#'
#' \code{travel} is a \code{bn.fit} object for the Bayesian network on a traveling preferences survey.
#'
#' @docType data
#'
#'
#' @usage travel
#'
#' @format The Bayesian network \code{travel} includes the following nodes:
#' \itemize{
#'   \item \bold{A}: three-level factor with levels \code{young}, \code{adult}, \code{old}. It indicates the age of an individual.
#'   \item \bold{S}: two level-factor with levels \code{M} (male) and \code{F} (female). It indicates the gender of an individual.
#'   \item \bold{E}: two level-factor with levels \code{high} and \code{uni}. It indicates the education level of an individual.
#'   \item \bold{O}: two level-factor with levels \code{emp} (employed) and \code{self} (self-employed). It indicates the occupation of an individual.
#'   \item \bold{R}: two level-factor with levels \code{small} and \code{big}. It indicates the size of the residence of an individual.
#'   \item \bold{T}: three level-factor with levels \code{car}, \code{train} and \code{other}. It indicates the preferred mean of transportation by an individual.
#' }
#'@source Scutari, M., & Denis, J. B. (2014). Bayesian networks: with examples in R. Chapman and Hall/CRC.
"travel"
