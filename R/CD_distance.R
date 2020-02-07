#' CD-distance
#'
#' \code{CD_distance} returns the Chan-Darwiche distance between a Bayesian network and its update after parameter variation.
#'
#' \cr
#' The Bayesian network on which parameter variation is being conducted should be expressed as a bn.fit object.
#' The name of the node to be varied, its level and its parent's level should be specified.
#' The parameter variation specified by the function is:\cr
#' \cr
#'  P ( \code{node} = \code{value_node} | parents = \code{value_parents} ) = \code{new_value} \cr
#'
#' @family dissimilarity measures
#'
#'@param bnfit object of class bn.fit
#'@param node character string. Node of which the conditional probability distribution is being changed.
#'@param value_node character string. Level of \code{node}.
#'@param value_parents character string. Levels of \code{node}'s parents. The levels should be defined according to the order of the parents in \code{bnfit[[node]][["parents"]]}. If \code{node} has no parents, then should be set to \code{NULL}.
#'@param new_value numeric vector with elements between 0 and 1. Values to which the parameter should be updated. It can take a specific value or more than one. In the case in which the user wants to consider more than one value,these should be defined through a vector with an increasing order of the elements. \code{new_value} can also take as value the character string \code{all}: in this case a sequence of possible parameter changes ranging from 0.05 to 0.95 will be considered.
#'@param covariation character string. Covariation scheme to be used for the updated Bayesian network. Can take values \code{uniform}, \code{proportional}, \code{orderp}, \code{all}.If equal to \code{all}, uniform, proportional and order-preserving co-variation schemes will be considered. Set by default to \code{proportional}.
#'@param plot boolean value. If \code{TRUE} the function returns a plot. If \code{covariation = "all"}, CD-distance for uniform (red), proportional (green), order-preserving (blue) co-variation schemes will be plotted.  Set by default to \code{TRUE}.
#'@param ... additional parameters to be added to the plot.

#'@importFrom stats coef
#'@importFrom graphics lines points
#'@importClassesFrom bnlearn bn.fit
#'@export
CD_distance <-
  function(bnfit,
           node,
           value_node,
           value_parents,
           new_value,
           covariation = "proportional",
           plot = TRUE,
           ...) {
    suppressWarnings(if (new_value == "all") {
      new_value2 <-
        sort(c(seq(0.05, 0.95, by = 0.05), coef(bnfit[[node]])[t(append(value_node, value_parents))]))
    } else{
      new_value2 <- new_value
    })
    new_value_op <- rep(NA, length(new_value2))
    bnfit.new <- vector("list", length = length(new_value2))
    if (covariation == "all") {
      bnfit.new.scheme <- vector("list", length = 3)
    } else{
      bnfit.new.scheme <- vector("list", length = 1)
    }
    if (covariation == "uniform" || covariation == "all") {
      for (i in 1:length(new_value2)) {
        bnfit.new[[i]] <- uniform_covar(
          bnfit = bnfit,
          node = node,
          value_node = value_node,
          value_parents = value_parents,
          new_value = new_value2[i]
        )
      }
      bnfit.new.scheme[[1]] <- bnfit.new

    }
    if (covariation == "proportional" || covariation == "all") {
      for (i in 1:length(new_value2)) {
        bnfit.new[[i]] <- proportional_covar(
          bnfit = bnfit,
          node = node,
          value_node = value_node,
          value_parents = value_parents,
          new_value = new_value2[i]
        )
      }
      if (covariation == "all") {
        bnfit.new.scheme[[2]] <- bnfit.new
      } else{
        bnfit.new.scheme[[1]] <- bnfit.new
      }
    }
    if (covariation == "orderp" || covariation == "all") {
      new_value_op <- new_value2
      length.node <- dim(coef(bnfit[[node]]))[1]
      node_values <- numeric(length.node)
      cpt.new <- coef(bnfit[[node]])
      for (i in 1:length.node) {
        node_values[i] <-
          cpt.new[t(append(unlist(dimnames(coef(
            bnfit[[node]]
          ))[1])[i], value_parents))]
      }
      order <- rank(node_values, ties.method = "first")
      s <-
        order[unname(which(unlist(dimnames(coef(
          bnfit[[node]]
        ))[1]) == value_node))]
      if (s == length.node & covariation == "orderp") {
        stop (
          "Order preserving covariation scheme can't work by varying the last parameter in the increasing order of the parameters."
        )
      }
      if (s == length.node & covariation == "all") {
        warning(
          "The CD distance with order-preserving scheme won't be returned since the last parameter in the increasing order of the parameters is being varied."
        )
        bnfit.new.scheme[[3]] <- NULL
      }
      if (s != length.node) {
        new_value_op <-
          new_value_op[new_value_op <= 1 / (1 + dim(coef(bnfit[[node]]))[1] - s)]
        if (length(new_value_op) != 0) {
          bnfit.new.op <- vector("list", length = length(new_value_op))
          for (i in 1:length(new_value_op)) {
            bnfit.new.op[[i]] <- orderp_covar(
              bnfit = bnfit,
              node = node,
              value_node = value_node,
              value_parents = value_parents,
              new_value = new_value_op[i]
            )
          }
          if (covariation == "all") {
            bnfit.new.scheme[[3]] <- bnfit.new.op
          } else{
            bnfit.new.scheme[[1]] <- bnfit.new.op
          }
        }
      }
    }
    if (covariation == "all") {
      if (length(bnfit.new.scheme) == 2) {
        CD <-
          data.frame(new_value2, rep(NA, length(new_value2)), rep(NA, length(new_value2)))
        colnames(CD) <- c("New value", "Uniform", "Proportional")
      } else{
        CD <-
          data.frame(new_value2,
                     rep(NA, length(new_value2)),
                     rep(NA, length(new_value2)),
                     rep(NA, length(new_value2)))
        colnames(CD) <-
          c("New value",
            "Uniform",
            "Proportional",
            "Order Preserving")
      }
    } else{
      CD <- data.frame(new_value2, rep(NA, length(new_value2)))
      colnames(CD) <- c("New value", "CD_distance")
    }

    if (length(new_value2) == length(new_value_op)) {
      for (k in 1:length(bnfit.new.scheme)) {
        for (j in 1:length(new_value2)) {
          quotients <- c()
          for (i in 1:dim(coef(bnfit[[node]]))[1]) {
            quotients[i] <-
              (coef(bnfit.new.scheme[[k]][[j]][[node]])[t(append(unlist(dimnames(
                coef(bnfit.new.scheme[[k]][[j]][[node]])
              )[1])[i], value_parents))]) / (coef(bnfit[[node]])[t(append(unlist(dimnames(
                coef(bnfit[[node]])
              )[1])[i], value_parents))])
          }
          CD[j, k + 1] <- log(max(quotients)) - log(min(quotients))
        }
      }
    } else{
      if (covariation == "all") {
        for (k in 1:2) {
          for (j in 1:length(new_value2)) {
            quotients <- c()
            for (i in 1:dim(coef(bnfit[[node]]))[1]) {
              quotients[i] <-
                (coef(bnfit.new.scheme[[k]][[j]][[node]])[t(append(unlist(dimnames(
                  coef(bnfit.new.scheme[[k]][[j]][[node]])
                )[1])[i], value_parents))]) / (coef(bnfit[[node]])[t(append(unlist(dimnames(
                  coef(bnfit[[node]])
                )[1])[i], value_parents))])
            }
            CD[j, k + 1] <-
              log(max(quotients)) - log(min(quotients))
          }
        }
        if (length(new_value_op) != 0) {
          for (t in 1:(length(new_value_op))) {
            quotients <- c()
            for (i in 1:dim(coef(bnfit[[node]]))[1]) {
              quotients[i] <-
                (coef(bnfit.new.scheme[[3]][[t]][[node]])[t(append(unlist(dimnames(
                  coef(bnfit.new.scheme[[3]][[t]][[node]])
                )[1])[i], value_parents))]) / (coef(bnfit[[node]])[t(append(unlist(dimnames(
                  coef(bnfit[[node]])
                )[1])[i], value_parents))])
            }
            CD[t, 4] <- log(max(quotients)) - log(min(quotients))
          }
        }
      } else{
        if (length(new_value_op) != 0) {
          for (t in 1:(length(new_value_op))) {
            quotients <- c()
            for (i in 1:dim(coef(bnfit[[node]]))[1]) {
              quotients[i] <-
                (coef(bnfit.new.scheme[[1]][[t]][[node]])[t(append(unlist(dimnames(
                  coef(bnfit.new.scheme[[1]][[t]][[node]])
                )[1])[i], value_parents))]) / (coef(bnfit[[node]])[t(append(unlist(dimnames(
                  coef(bnfit[[node]])
                )[1])[i], value_parents))])
            }
            CD[t, 2] <- log(max(quotients)) - log(min(quotients))
          }
        }
      }
    }
    if (all(is.na(CD[,-1]))) {
      plot <- FALSE
      warning("The plot won't be showed since all the values are not possible")
    }
    if (plot == TRUE) {
      if (covariation == "all") {
        if (nrow(CD) == 1) {
          plot <-
            plot(
              x = CD[, 1],
              y = CD[, 2],
              main = "CD-distance",
              xlab = "New_value",
              ylab = "CD-distance",
              col = "red",
              pch = 16,
              xlim = c(0, 1),
              ylim = c(min(unlist(CD[, -1])[is.finite(unlist(CD[, -1]))]), max(unlist(CD[, -1])[is.finite(unlist(CD[, -1]))])),
              ...
            )
          points(
            x = CD[, 1],
            y = CD[, 3],
            col = "green",
            pch = 16,
            ...
          )
          if (ncol(CD) == 4) {
            points(
              x = CD[, 1],
              y = CD[, 4],
              col = "blue",
              pch = 16,
              ...
            )
          }
        } else{
          plot <-
            plot(
              x = CD[, 1],
              y = CD[, 2],
              main = "CD-distance",
              xlab = "New_value",
              ylab = "CD-distance",
              col = "red",
              xlim = c(0, 1),
              ylim = c(min(unlist(CD[, -1])[is.finite(unlist(CD[, -1]))]), max(unlist(CD[, -1])[is.finite(unlist(CD[, -1]))])),
              type = "l",
              ...
            )
          lines(x = CD[, 1],
                y = CD[, 3],
                col = "green",
                ...)
          if (ncol(CD) == 4) {
            lines(x = CD[, 1],
                  y = CD[, 4],
                  col = "blue",
                  ...)
          }
        }
      } else{
        if (nrow(CD) == 1) {
          plot <-
            plot(
              x = CD[, 1],
              y = CD[, 2],
              main = "CD-distance",
              xlab = "New_value",
              ylab = "CD-distance",
              xlim = c(0, 1),
              ylim = c(min(unlist(CD[, -1])[is.finite(unlist(CD[, -1]))]), max(unlist(CD[, -1])[is.finite(unlist(CD[, -1]))])),
              pch = 16,
              ...
            )
        } else{
          plot <-
            plot(
              x = CD[, 1],
              y = CD[, 2],
              main = "CD-distance",
              xlab = "New_value",
              ylab = "CD-distance",
              type = "l",
              xlim = c(0, 1),
              ylim = c(min(unlist(CD[, -1])[is.finite(unlist(CD[, -1]))]), max(unlist(CD[, -1])[is.finite(unlist(CD[, -1]))])),
              ...
            )
        }
      }
    }

    return(list(CD = CD, plot = plot))
  }
