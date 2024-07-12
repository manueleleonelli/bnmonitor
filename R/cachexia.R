#' Bayesian networks for a cachexia study
#'
#' Continuous Bayesian networks comparing the dependence of metabolomics for people who suffer and do not suffer of Cachexia
#'
#' @docType data
#'
#'
#' @format Continuous Bayesian networks over six metabolomics: Adipate (A), Betaine (B), Fumarate (F), Glucose (GC), Glutamine (GM) and Valine (V).
#'     The networks \code{cachexia_gbn} and \code{cachexia_ci} are for people suffering of cachexia and of class \code{GBN} and \code{CI} respectively.
#'     The networks \code{control_gbn} and \code{control_ci} are for people not suffering of cachexia and of class \code{GBN} and \code{CI} respectively.
#'     The original dataset is stored in \code{cachexia_data}.
#' @aliases cachexia_gbn
#' @name cachexia
#'
#'@references	C. Görgen & M. Leonelli (2020), Model-preserving sensitivity analysis for families of Gaussian distributions.  Journal of Machine Learning Research, 21: 1-32.
#'
NULL

#'
#' @rdname cachexia
"cachexia_gbn"

#' @rdname cachexia
"cachexia_ci"

#'
#' @rdname cachexia
"control_gbn"

#' @rdname cachexia
"control_ci"

#' @rdname cachexia
"cachexia_data"

