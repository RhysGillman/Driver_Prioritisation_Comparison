
#' @title Calculating driving force matrix
#' @description
#'   Calculating driving force matrix for each patient
#' @param network the edge list of signal transduction network.
#' @param diffexpgene a list consists with P-value (or P-adjust) and log2 fold change matrix
#' which indicates whether a gene is differentially expressed.
#' @param enrichreactomeRes abnormal pathways for each patient.
#' @param reactome Reactome pathway information.
#' @param parallelworker number of workers that parallel calculating differentally
#' expression genes of every patients.
#' @param th.path giving threshold of enrichment pathways.
#' @param th.FC giving threshold of differentially expression log2 fold change.
#' @param th.enrich.type indicating whether P-value or P-adjust used to calculate
#' enrichment p-value of upstream regulators
#' @param th.enrich giving threshold of enrichment p-value.
#' @param th.zscore giving threshold of causality z-score.
#' @importFrom clusterProfiler enricher
#' @import doParallel
#' @import foreach
#' @import parallel
#'
#' @return A list consists two elements, one of them consists with driving force matrix of every patients,
#' another consists with individual z-score.
#' @export calDrivingForce
#'
#' @examples
#' \dontrun{
#'drivingforcematrix  <- calDrivingForce(network = STN,
#'                                       diffexpgene = diffexpgene,
#'                                       enrichreactomeRes = enrichreactomeRes,
#'                                       reactome = reactome,
#'                                       parallelworker = 15)
#' }
#' @author Yan Li

calDrivingForce <- function(network,
                            diffexpgene,
                            enrichreactomeRes,
                            reactome,
                            parallelworker,
                            th.path = 0.05,
                            th.FC = 3,
                            th.enrich.type = 'p.adj',
                            th.enrich = 0.01,
                            th.zscore = 2){

  network <- rbind(network, data.frame(Source_nodes=unique(c(network$Source_nodes,network$Target_nodes)),
                                       Target_nodes=unique(c(network$Source_nodes,network$Target_nodes)),
                                       Score = rep(0, length(unique(c(network$Source_nodes,network$Target_nodes)))),
                                       Source = rep('self-active',
                                                    length(unique(c(network$Source_nodes,network$Target_nodes))))))
  network$Score[network$Source == 'Others'] <-
    network$Score[network$Source == 'Others'] * 0.5
  nodename <- unique(c(network$Source_nodes, network$Target_nodes))

  # Constructing adjacent matrix of STN aij(i -> j)
  signalNetadj <- matrix(0, nrow = length(nodename), ncol = length(nodename),
                         dimnames = list(nodename, nodename))
  signalNetadj[match(network$Source_nodes, nodename) +
                 (match(network$Target_nodes, nodename) - 1) * length(nodename)] <- network$Score
  signalNetadjnorm <- signalNetadj %*% diag(1 / (colSums(abs(signalNetadj)) + 0.0001))

  gene.name <- rownames(diffexpgene$diffFC)
  reactomesplit <- split(reactome, reactome$V4)

  if(parallelworker > 1) {
    cl <- makeCluster(parallelworker)
    registerDoParallel(cl)
  }
  genescore_individual <- foreach (i = 1:dim(diffexpgene$diffFC)[2]) %dopar% {
    enrichgene <- strsplit(enrichreactomeRes[[i]]$geneID[enrichreactomeRes[[i]]$pvalue <= th.path], '/')
    enrichname <- enrichreactomeRes[[i]]$ID[enrichreactomeRes[[i]]$pvalue <= th.path]
    pathpvalue <- enrichreactomeRes[[i]]$pvalue[enrichreactomeRes[[i]]$pvalue <= th.path]
    names(enrichgene) <- enrichname

    # Calculating z-score for upstream regulators of abnormal pathway DEGs
    res <- list()
    for (j in 1:length(enrichname)) {
      if (is.na(th.FC)) {
        flag <- (gene.name %in% reactomesplit[[which(names(reactomesplit) == enrichname[j])]]$V1) |
          (gene.name %in% enrichgene[[j]])
      }
      else {
        flag <- abs(diffexpgene$diffFC[, i]) >= th.FC &
          (gene.name %in% reactomesplit[[which(names(reactomesplit) == enrichname[j])]]$V1) |
          (gene.name %in% enrichgene[[j]])
      }
      if (sum(flag) < 1) next
      IG <- gene.name[flag]
      if (sum(IG %in% intersect(gene.name, unique(c(network$Source_nodes,network$Target_nodes)))) == 0) next
      tmp <- clusterProfiler::enricher(IG, TERM2GENE = network[, c(1, 2)],
                                       universe = intersect(gene.name, unique(c(network$Source_nodes,
                                                                                network$Target_nodes))),
                                       minGSSize = 1, maxGSSize = 1500)@result
      fcdir <- sign(diffexpgene$diffFC[, i][match(IG, gene.name)])
      zscore <- rep(0, dim(tmp)[1])
      for (k in 1:dim(tmp)[1]) {
        currentedgelist <- network[(network$Source_nodes %in% tmp[k,1] & network$Target_nodes %in% IG), ]
        if ('0' %in% currentedgelist$Score) {
          overlapscore <- sum(sign(currentedgelist$Score) * fcdir[match(currentedgelist$Target_nodes, IG)])
          zscore[k] <- (overlapscore + sign(overlapscore)) / sqrt(dim(currentedgelist)[1])
        } else {
          zscore[k] <- sum(sign(currentedgelist$Score) * fcdir[match(currentedgelist$Target_nodes, IG)]) /
            sqrt(dim(currentedgelist)[1])}}
      tmp$pathpvalue <- rep(pathpvalue[j], length(zscore))
      tmp$zscore <- zscore
      res[[j]] <- tmp
    }
    names(res) <- enrichname[1:length(res)]
    res <- res[!(sapply(res, is.null))]
    res <- lapply(res, function(x){
      if (th.enrich.type == 'p.adj') enrich.p <- x[,6]
      else if (th.enrich.type == 'p.value') enrich.p <- x[,5]
      res <- x[enrich.p < th.enrich & abs(x[, 11]) > th.zscore & !(is.nan(x[,11])),]
      return(res)
    })
    res <- res[sapply(res, function(x) return(dim(x)[1])) != 0]
    if (length(res) == 0) return(NULL)
    enrichgene <- enrichgene[match(names(res), names(enrichgene))]

    # Network attribution
    genescore <- matrix(0, dim(signalNetadj)[1], length(res), dimnames = list(nodename, names(res)))
    for (j in 1:length(res)) {
      if (is.na(th.FC)) {
        flag <- (gene.name %in% reactomesplit[[which(names(reactomesplit) == names(res)[j])]]$V1) |
          (gene.name %in% enrichgene[[j]])
      }
      else {
        flag <- abs(diffexpgene$diffFC[, i]) >= th.FC &
          (gene.name %in% reactomesplit[[which(names(reactomesplit) == names(res)[j])]]$V1) |
          (gene.name %in% enrichgene[[j]])
      }

      IG <- gene.name[flag]
      IG <- IG[IG %in% nodename]
      fcdir <- sign(diffexpgene$diffFC[, i][match(IG, gene.name)])
      subnet <- signalNetadj[match(res[[j]]$ID, rownames(signalNetadj)), match(IG, colnames(signalNetadj))]
      if (is.null(dim(subnet))) subnet <- subnet / (abs(subnet) + 0.0001)
      else subnet <- subnet %*% diag(1 / (colSums(abs(subnet)) + 0.0001))
      inits <- subnet %*% matrix(fcdir / length(fcdir), length(fcdir), 1)
      genescore[match(res[[j]]$ID, rownames(genescore)),j] <- inits
    }
    genescore <- genescore + signalNetadjnorm %*% genescore + signalNetadjnorm %*% signalNetadjnorm %*% genescore +
      signalNetadjnorm %*% signalNetadjnorm %*% signalNetadjnorm %*% genescore
    return(list(genescore = genescore, zscore_indiv = res))
  }
  if(parallelworker > 1) stopCluster(cl)
  genescore <- lapply(genescore_individual, function(x){
    if (is.null(x)) return(x)
    else return(x[[1]])
  })
  zscore <- lapply(genescore_individual, function(x){
    if (is.null(x)) return(x)
    else return(x[[2]])
  })
  names(genescore) <- colnames(diffexpgene$diffFC)
  names(zscore) <- colnames(diffexpgene$diffFC)
  return(list(genescore = genescore, zscore = zscore))
}
