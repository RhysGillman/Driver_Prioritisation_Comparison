
#' @title Get driver genes
#' @description
#'   Get individual-level and population-level driver genes from driving force matrix.
#' @param drivingforcelist A list consists with driving force matrix of each patient.
#' @param cancerData A list consists with expression and mutation data after per-processing.
#' @param topgenes_ind Number of top genes that considered as individual driver genes.
#'
#' @return A list with individual driver list and ranked population driver list
#' @export getDriverGenes
#'
#' @examples
#' \dontrun{
#' drivergenes <- getDriverGenes(drivingforcelist = genescore_individual,
#'                               cancerData = cancerData)
#' }
#' @author Yan Li

getDriverGenes <- function(drivingforcelist,
                           cancerData,
                           topgenes_ind = 5){
  maf <- cancerData$mutationData
  genescore_individual <- drivingforcelist[!(sapply(drivingforcelist, is.null))]
  sample_snp  <- sapply(strsplit(maf$Tumor_Sample_Barcode, split = '-'),
                        function(x) return(paste(x[1:4], collapse = '-')))
  sample_score <- sapply(strsplit(names(genescore_individual), split = '-'),
                         function(x) return(paste(x[1:4], collapse = '-')))
  genescore_individual_snp <- genescore_individual
  for (i in 1:length(genescore_individual_snp)) {
    if (sample_score[i] %in% sample_snp) {
      temp <- maf[sample_snp %in% sample_score[i],]
      genescore_individual_snp[[i]] <- genescore_individual_snp[[i]][rownames(genescore_individual_snp[[i]]) %in%
                                                                       temp$Hugo_Symbol,]
    }
    else genescore_individual_snp[[i]] <- 1
  }
  genescore_individual_snp <- genescore_individual_snp[sample_score %in% sample_snp]
  for (i in 1:length(genescore_individual_snp)) {
    if (is.null(dim(genescore_individual_snp[[i]]))) next
    else genescore_individual_snp[[i]] <- rowSums(genescore_individual_snp[[i]])
  }
  genescore_individual_snp <- lapply(genescore_individual_snp,
                                     function(x) return(x[order(abs(x),decreasing = T)]))

  drivergene_ind <- lapply(genescore_individual_snp, function(x) return(x[1:topgenes_ind]))
  topgenes <- drivergene_ind
  names(topgenes) <- NULL
  drivergene_pop <- as.data.frame(sort(table(names(unlist(topgenes))), decreasing = T))
  return(list(drivergene_ind = drivergene_ind, drivergene_pop = drivergene_pop))
}
