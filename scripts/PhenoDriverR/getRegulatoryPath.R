
#' @title Get regulatory paths from driver gene to abnormal pathways
#' @description
#'   Reconstructing sub-network which describes the possible regulatory paths from interest driver genes
#'   to its associated abnormal pathways.
#' @param drivergenes A list with individual driver list and ranked population driver list.
#' @param zscore_ind A list consist with individual z-score.
#' @param network the edge list of signal transduction network.
#' @param quary indicating the query driver gene
#' @param quary_pathway indicating which abnormal pathways need to be constructed regulatory paths
#' @import igraph
#'
#' @return edge list of regulatory paths from query driver gene to query abnormal pathways, which can be
#' used directly as the Cytoscape input to visualize this sub-network.
#' @export getRegulatoryPath
#'
#' @examples
#' \dontrun{
#' regulation_path_edgelist <- getRegulatoryPath(drivergenes = drivergenes,
#'                                               zscore_ind = genescore_individual$zscore,
#'                                               network = STN,
#'                                               quary = 'TP53',
#'                                               quary_pathway = c('EML4 and NUDC in mitotic spindle formation',
#'                                                                 'Resolution of Sister Chromatid Cohesion',
#'                                                                 'Mitotic Prometaphase',
#'                                                                 'RHO GTPases Activate Formins',
#'                                                                 'Separation of Sister Chromatids',
#'                                                                 'Amplification of signal from unattached kinetochores via a MAD2 inhibitory signal'))
#' }
#' @author Yan Li

getRegulatoryPath <- function(drivergenes,
                              zscore_ind,
                              network,
                              quary,
                              quary_pathway){
  network[network$source %in% c('Reactom', 'complex'), 3] <-
    network[network$source %in% c('Reactom', 'complex'), 3] * 0.5
  nodename <- unique(c(network$Source_nodes, network$Target_nodes))
  signalNetadj <- matrix(0, nrow = length(nodename), ncol = length(nodename), dimnames = list(nodename, nodename))
  signalNetadj[match(network$Source_nodes, nodename) + (match(network$Target_nodes, nodename) - 1) * length(nodename)] <- network$Score
  signalNetadj <- abs(signalNetadj)
  signalNet <- graph_from_adjacency_matrix(signalNetadj, 'directed', weighted = T)

  s <- sapply(drivergenes$drivergene_ind, function(x) return(quary %in% names(x)))
  quary_samples <- names(s)[s == T]
  zscore_ind <- zscore_ind[names(zscore_ind) %in% quary_samples]

  regulation_path_indiv <- lapply(zscore_ind, .get_regulation_path, quary_pathway, signalNet, quary)
  regulation_path_merge <- do.call(rbind, regulation_path_indiv)
  regulation_path_merge$lable <- paste(regulation_path_merge$source, regulation_path_merge$target, sep = '|')
  regulation_path_merge$score <- network$Score[match(regulation_path_merge$lable, paste(network$Source_nodes, network$Target_nodes, sep = '|'))]
  regulation_path_merge$db <- network$Source[match(regulation_path_merge$lable, paste(network$Source_nodes, network$Target_nodes, sep = '|'))]
  regulation_path_edgelist <- unique(regulation_path_merge)
  regulation_path_edgelist$freq <- table(regulation_path_merge$lable)[match(regulation_path_edgelist$lable,
                                                                            names(table(regulation_path_merge$lable)))]
  regulation_path_edgelist <- regulation_path_edgelist[!(is.na(regulation_path_edgelist$source)) & !(is.na(regulation_path_edgelist$target)), ]
  return(regulation_path_edgelist)
}


.convert_path_2_df <- function(path) {
  path <- as_ids(path)
  if(length(path) > 4 | length(path) == 1) return(NA)
  else return(data.frame(source = path[1:(length(path) - 1)], target = path[2:length(path)]))
}

.re_target <- function(path) {
  path <- as_ids(path)
  if(length(path) > 4 | length(path) == 1) return(NA)
  else return(path[length(path)])
}

.get_regulation_path <- function(res, quary_pathway, signalNet, quary){
  res <- res[names(res) %in% quary_pathway]
  if (length(res) == 0) return(NA)
  for (i in 1:length(res)) {
    res[[i]]$pathname <- rep(names(res)[i], dim(res[[i]])[1])
  }
  path_df <- lapply(res, function(x, signalNet) {
    netpath <- all_shortest_paths(signalNet, quary,
                                  x$ID[x$Count / as.numeric(sapply(strsplit(x$GeneRatio, '/'),
                                                                   function(x) return(x[2]))) > 0.3], 'out', NA)
    tmp <- lapply(netpath$res, .convert_path_2_df)
    tmp <- unique(do.call(rbind, tmp))
    targetnode <- unique(sapply(netpath$res, .re_target))
    tmp <- rbind(tmp, data.frame(source = targetnode, target = rep(unique(x$pathname),length(targetnode))))
    return(tmp)
  }, signalNet)
  names(path_df) <- NULL
  path_df <- unique(do.call(rbind, path_df))
  return(path_df)
}
