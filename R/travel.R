#' Bayesian network on travel survey
#'
#' \code{travel} is a bn.fit object for the Bayesian network on a travelling preferences survey.
#'
#' @docType data
#'
#' @family datasets
#'
#' @usage travel
#'
#' @format The Bayesian network \code{travel} comprehends the following nodes:
#' \itemize{
#'   \item \bold{A}: three-level factor with levels \emph{young}, \emph{adult}, \emph{old}. Indicates the age of a person.
#'   \item \bold{S}: two level-factor with levels \emph{M} (male) and \emph{F} (female). Indicates the sex of a person.
#'   \item \bold{E}: two level-factor with levels \emph{high} and \emph{uni}. Indicates the education level of a person.
#'   \item \bold{O}: two level-factor with levels \emph{emp} (employed) and \emph{self} (self-employed). Indicates the occupation of a person.
#'   \item \bold{R}: two level-factor with levels \emph{small} and \emph{big}. Indicates the size of the residence of a person.
#'   \item \bold{T}: three level-factor with levels \emph{car}, \emph{train} and \emph{other}. Indicates the preferred mean of transport by a person.
#' }
#'@source Marco Scutari, Jean-Baptiste Denis (2014). "Bayesian Networks:With Examples in R".
"travel"
