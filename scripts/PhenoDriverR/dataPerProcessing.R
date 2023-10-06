
#' @title Data per-processing for gene expression data and mutation data downloaded from TCGA
#' @description
#'   Per-processing for raw counts gene expression data and mutation data,
#'   and filter out samples without matched gene expression data and mutation data
#' @param expressionData a list contains gene expression matrix for NT and TP samples,
#' which named as 'Normal' and 'Tumor' respactively.
#' @param mutationData a dataframe contains mutation terms for each patients.
#' @param annotationcol indicate the number of column of gene annotations.
#' @param removeDuplicates indicate the column of gene annotation that used to
#' remove duplicate genes. If there is no need to perform de-duplication, set
#' removeDuplicates=NULL.
#' @param removeSilenceSNP Should remove silence SNPs from mutation data?
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' @importFrom stats median
#'
#' @return A list with the results after per-processing
#' @export dataPerProcessing
#'
#' @examples
#' \dontrun{
#' cancerData <- dataPerProcessing(expressionData = expressionData,
#'                                 mutationData = mutationData,
#'                                 annotationcol = 3,
#'                                 removeDuplicates = 2,
#'                                 removeSilenceSNP = F)
#' }
#' @author Yan Li

dataPerProcessing <- function(expressionData,
                              mutationData,
                              annotationcol,
                              removeDuplicates = NULL,
                              removeSilenceSNP = F) {

  allexp <- merge(expressionData$Normal, expressionData$Tumor)
  col.factor <- estimateSizeFactorsForMatrix(allexp[, -1:(-1*annotationcol)])
  allexp[, -1:(-1*annotationcol)] <- as.data.frame(t(t(allexp[, -1:(-1*annotationcol)]) * (1 / col.factor)))

  # Removing duplicate genes
  if (!(is.null(removeDuplicates))) {
    gene.name <- allexp[, removeDuplicates]
    repetition <- as.data.frame(table(gene.name))
    repetition <- repetition[repetition$Freq > 1, ]
    del <- c()
    for (i in 1:dim(repetition)[1]) {
      rownum <- which(gene.name == repetition[i, 1])
      allexp[rownum[1], -1:(-1*annotationcol)] <-
        apply(as.matrix(allexp[rownum, -1:(-1*annotationcol)]), 2, median)
      del <- c(del, rownum[2:length(rownum)] * -1)
    }
    allexp <- allexp[del, ]
    gene.name <- gene.name[del]
  }

  # Removing genes with low expression values
  NTindex <- match(colnames(expressionData$Normal)[-1:(-1*annotationcol)],
                   colnames(allexp))
  TPindex <- match(colnames(expressionData$Tumor)[-1:(-1*annotationcol)],
                   colnames(allexp))
  allexp <- allexp[(rowSums(allexp[NTindex] > 20) >= (0.8 * length(NTindex))) |
                     (rowSums(allexp[TPindex] > 20) >= (0.8 * length(TPindex))), ]
  expressionData$Normal <- allexp[, c(1:annotationcol, NTindex)]
  expressionData$Tumor <- allexp[, c(1:annotationcol, TPindex)]

  # Removing hyper-mutated patients
  mutationData <-
    mutationData[mutationData$Tumor_Sample_Barcode %in%
                              names(table(mutationData$Tumor_Sample_Barcode)
                                    [table(mutationData$Tumor_Sample_Barcode) <= 1000]), ]

  # Removing silence mutations
  if (removeSilenceSNP) {
    mutationData <- mutationData[
      mutationData$Variant_Classification != 'Silent', ]
  }

  cancerData <- list(expressionData = expressionData,
                     mutationData = mutationData)
  return(cancerData)
}
