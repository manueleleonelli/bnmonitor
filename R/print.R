
#' Printing methods
#'
#' Printing methods for outputs of \code{bnmonitor} functions
#'
#' @param x an appropriate object
#' @param ... for compatibility
#' @name print
#' @return Printing specific to the object it is applied to.
NULL


#'@export
#'@rdname print
#'
#'
print.sensitivity <- function(x,...){
  print(x$sensitivity)
  invisible(x)
}

#'@export
#'@rdname print
#'
#'
print.diameter <- function(x,...){
  print(x$Diameter)
  invisible(x)
}

#'@export
#'@rdname print
#'
#'
print.mutualinfo <- function(x,...){
  print(x$MutualInfo)
  invisible(x)
}

#'@export
#'@rdname print
#'
#'
print.dwi <- function(x,...){
  print(x$DWI)
  invisible(x)
}

#'@export
#'@rdname print
#'
#'
print.ewi <- function(x,...){
  print(x$EWI)
  invisible(x)
}

#'@export
#' @rdname print
print.kl <- function(x,...){
  print(x$KL)
  invisible(x)
}


#'@export
#'
#'@rdname print
#'
print.CD <- function(x,...){
  print(x$CD)
  invisible(x)
}


#'@export
#'
#'
#'@rdname print
#'
print.fro <- function(x,...){
  print(x$Frobenius)
  invisible(x)
}


#' @importClassesFrom bnlearn bn.fit
#'@export
#'@rdname print
#'
print.node_monitor <- function(x,...){
  print(x$Global_Monitor)
  invisible(x)
}


#'@export
#' @rdname print
#'
print.jeffreys <- function(x,...){
  print(x$Jeffreys)
  invisible(x)
}


#' @importClassesFrom bnlearn bn.fit
#'@export
#'@rdname print
#'
print.final_node_monitor <- function(x,...){
  print(x$Node_Monitor)
  invisible(x)
}



#'@export
#'@rdname print
#'
print.seq_cond_monitor <- function(x,...){
  temp <- x$Seq_Cond_Monitor
  temp <- temp[is.finite(temp)]
  cat("Conditional Node Monitor for", x$node.name,"\n",
      "Minimum ", "\t", min(temp,na.rm = TRUE), "\n",
      "Maximum", "\t", max(temp,na.rm = TRUE))
  invisible(x)
}


#'@export
#'
#'@rdname print
#'
print.seq_pa_ch_monitor <- function(x,...){
  temp <- x[1:length(x)]
  temp <- temp[is.finite(temp)]
  cat("Parent Child Node Monitor","\n",
      "Minimum ", "\t", min(temp,na.rm = TRUE), "\n",
      "Maximum", "\t", max(temp,na.rm = TRUE))
  invisible(x)
}

#'@export
#'
#'@rdname print
print.seq_marg_monitor <- function(x,...){
  temp <- x$Seq_Marg_Monitor
  temp <- temp[is.finite(temp)]
  cat("Marginal Node Monitor for", x$node.name,"\n",
      "Minimum ", "\t", min(temp,na.rm = TRUE), "\n",
      "Maximum", "\t", max(temp,na.rm = TRUE))
  invisible(x)
}

