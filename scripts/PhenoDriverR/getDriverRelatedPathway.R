
#' @title Get driver related pathway/phenotype
#' @description
#'   Calculating driver-associated abnormal pathways for a giving driver gene
#' @param drivingforcelist A list consists with driving force matrix of each patient.
#' @param enrichreactomeRes Abnormal pathways for each patient.
#' @param drivergenes A list with individual driver list and ranked population driver list.
#' @param quary indicating the query driver gene
#'
#' @return A dataframe contains the information of driver-associated pathways
#' @export getDriverRelatedPathway
#'
#' @examples
#' \dontrun{
#' tp53relatedpathway <- getDriverRelatedPathway(drivingforcelist = genescore_individual,
#'                                               enrichreactomeRes = enrichreactomeRes,
#'                                               drivergenes = drivergenes,
#'                                               quary = 'TP53')
#' }
#' @author Yan Li

getDriverRelatedPathway <- function(drivingforcelist,
                                    enrichreactomeRes,
                                    drivergenes,
                                    quary){
  genescore_individual <- drivingforcelist[!(sapply(drivingforcelist, is.null))]
  genescore_single <- lapply(genescore_individual, function(x) return(x[rownames(x) == quary,]))
  path_name <- unique(unlist(lapply(drivingforcelist, colnames)))
  genescore_matrix <- matrix(0, length(genescore_single), length(path_name), dimnames = list(names(genescore_single), path_name))
  for (i in 1:length(genescore_single)) {
    genescore_matrix[i, match(names(genescore_single[[i]]), path_name)] <- genescore_single[[i]]
  }
  s <- sapply(drivergenes$drivergene_ind, function(x) return(quary %in% names(x)))
  s <- s[s == TRUE]

  term2index <- apply(genescore_matrix, 2, function(x, samplename){
    return(samplename[x != 0])
  }, samplename = rownames(genescore_matrix))
  term2index <- term2index[sapply(term2index, length)!=0]
  term <- rep(names(term2index), sapply(term2index, length))
  names(term2index) <- NULL
  index <- unlist(term2index)
  term2index <- data.frame(term = term, index = index)
  driverrelatedpathway <- clusterProfiler::enricher(names(s), TERM2GENE = term2index,
                                                    universe = rownames(genescore_matrix),
                                                    minGSSize = 5, maxGSSize = 1500)@result
  pathway_NP <- colSums(genescore_matrix)
  pathway_NP[pathway_NP > 0] <- 'UP'
  pathway_NP[pathway_NP < 0] <- 'DOWN'
  driverrelatedpathway$`UP/DOWN` <- pathway_NP[match(driverrelatedpathway$ID, names(pathway_NP))]
  return(driverrelatedpathway)
}
