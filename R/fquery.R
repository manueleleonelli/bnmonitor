#'@import bnlearn
#'@importClassesFrom bnlearn bn.fit
#'@importFrom stats coef
#'@importFrom bnlearn as.grain
#'@importFrom gRain setEvidence
#'@importFrom gRain querygrain
#'@importFrom gRbase compile

fquery <-
  function (x,
            new_value,
            bnfit,
            interest_node,
            interest_node_value,
            evidence_nodes = NULL,
            evidence_states = NULL,
            node,
            value_node,
            value_parents) {
    bnfit.new <- tryCatch(
      proportional_covar(
        bnfit = bnfit,
        node = node,
        value_node = value_node,
        value_parents = value_parents,
        new_value = x
      )
    )
    grain.bn <- compile(as.grain(bnfit.new))
    suppressWarnings(if (!is.null(evidence_nodes)) {
      grain.bn <- setEvidence(grain.bn, nodes = evidence_nodes,
                              states = evidence_states)
    })
    fquery <-
      querygrain(grain.bn, nodes = interest_node)[[interest_node]][interest_node_value] -
      new_value
    return(fquery)
  }
