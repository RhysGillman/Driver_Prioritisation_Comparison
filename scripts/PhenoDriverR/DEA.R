
#' @title Personalized differential expression analysis (DEA) based on noise model
#' @description
#'   Differential exprression analysis (DEA) for each patient
#' @param exp a list consists with two elements, expression matrix of NT samples
#' and TP samples. Each elements is a data frame with gene as rows and patient as
#' colnums.
#' @param parallelworker number of workers that parallel calculating differentally
#' expression genes of every patients.
#' @param genenamecol indicate the anotation column that regarded as gene names.
#' @param annotationcol indicate the number of column of gene annotations.
#' @importFrom gtools combinations
#' @importFrom stats p.adjust
#' @import doParallel
#' @import foreach
#' @import parallel
#'
#' @return A list consists with the differential expression p-value, p-adj and
#' log2 fold change.
#' @export DEA
#'
#' @examples
#' \dontrun{
#' diffexpgene <- DEA(exp = cancerData$expressionData,
#'                    parallelworker = 4,
#'                    genenamecol = 2,
#'                    annotationcol = 2)
#' }
#' @author Yan Li

cal.pvalue <- function(x, numeric.model) {
  distrubution <- numeric.model[[as.character(x[1])]]
  if (is.null(distrubution)) return(c(1, sign(x[2] - log(x[1] + 1))))
  if (x[2] >= log(x[1] + 1)) return(c(sum(distrubution$V2 >= x[2]) / dim(distrubution)[1], 1))
  else return(c(sum(distrubution$V2 <= x[2]) / dim(distrubution)[1], -1))
}


DEA <- function(exp,
                parallelworker,
                genenamecol = NA,
                annotationcol = NA){

  if(!(is.na(genenamecol))) {
    rownames(exp$Normal) <- exp$Normal[, genenamecol]
    rownames(exp$Tumor) <- exp$Tumor[, genenamecol]
  }
  if(!(is.na(annotationcol))) {
    exp$Normal <- exp$Normal[, -1:(-1*annotationcol)]
    exp$Tumor <- exp$Tumor[, -1:(-1*annotationcol)]
  }

  comb <- combinations(dim(exp$Normal)[2], 2)
  p <- matrix(0, dim(comb)[1] * dim(exp$Normal)[1], 2)
  for (i in 1:dim(comb)[1]) {
    p[((i-1)*dim(exp$Normal)[1]+1):(i*dim(exp$Normal)[1]),] <-
      as.matrix(exp$Normal[,comb[i,]])
  }

  p <- rbind(p,p[,c(2,1)])
  p <- round(p*10)
  p[,2] <- log(p[,2] + 1)
  a <- gc()
  numeric.model <- split(as.data.frame(p), as.factor(p[,1]))
  p <- 0
  a <- gc()

  if(parallelworker > 1) {
    cl <- makeCluster(parallelworker)
    registerDoParallel(cl)
  }
  diffRes <- foreach (i = 1:dim(exp$Tumor)[2], .export = c("cal.pvalue")) %dopar% {
    test <- data.frame(x=round(apply(exp$Normal, 1, median)*10),y=log(exp$Tumor[, i]*10 + 1))
    res <- apply(test, 1, cal.pvalue, numeric.model)
    pvalue <- res[1,]
    sign <- res[2,]
    padj <- p.adjust(pvalue, 'BH')
    return(list(pvalue = pvalue, sign = sign, padj = padj))
  }
  if(parallelworker > 1) stopCluster(cl)

  names(diffRes) <- colnames(exp$Tumor)
  diffexpgene <- list()
  diffexpgene$diffpadj <- sapply(diffRes, function(x) return(x$padj))
  diffexpgene$diffpval <- sapply(diffRes, function(x) return(x$pvalue))
  diffexpgene$diffsign <- sapply(diffRes, function(x) return(x$sign))
  diffexpgene$diffFC <- apply(exp$Tumor,2,function(x,y) return(log2((x + 0.1) / (y + 0.1))),
                              apply(exp$Normal, 1, median))
  return(diffexpgene)
}


