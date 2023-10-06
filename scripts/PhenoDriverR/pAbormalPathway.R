
#' @title Abnormal pathways discovering
#' @description
#'   Discovering individual patient's clinically abnormal phenotypes (pathways)
#' @param diff a list consists with P-value (or P-adjust) matrix which indicates whether
#' a gene is differentially expressed.
#' @param reactome Reactome pathway information
#' @param parallelworker number of workers that parallel calculating differentally
#' expression genes of every patients.
#' @param th.diff.type indicating whether P-value or P-adjust used to ditect DEGs.
#' @param th.diff giving threshold of differentially expression pvalue.
#' @importFrom clusterProfiler enricher
#' @importFrom stats p.adjust
#' @import doParallel
#' @import foreach
#' @import parallel
#'
#' @return  the results of enrichment analysis for every patient.
#' @export pAbnormalPathway
#'
#' @examples
#' \dontrun{
#' abnormalPathway <- pAbnormalPathway(diff = diffexpgene,
#'                                     reactome = reactome,
#'                                     parallelworker = 4)
#' }
#' @author Yan Li

pAbnormalPathway <- function(diff,
                             reactome,
                             parallelworker,
                             th.diff.type = 'p.value',
                             th.diff = 0.01){
  if (th.diff.type == 'p.value') allpvalue <- diff$diffpval
  else if (th.diff.type == 'p.adj') allpvalue <- diff$diffpadj
  if(parallelworker > 1) {
    cl <- makeCluster(parallelworker)
    registerDoParallel(cl)
  }
  enrichreactomeRes <- foreach(i = 1:dim(allpvalue)[2])  %dopar% {
    tmp <- clusterProfiler::enricher(rownames(allpvalue)[abs(allpvalue[, i]) < th.diff],
                                     TERM2GENE = reactome[, c(4, 1)],
                                     universe = rownames(allpvalue),
                                     minGSSize = 5,
                                     maxGSSize = 1500)
    return(tmp@result)
  }
  if(parallelworker > 1) stopCluster(cl)
  names(enrichreactomeRes) <- colnames(allpvalue)

  return(enrichreactomeRes)
}
