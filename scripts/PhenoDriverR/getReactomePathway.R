
#' @title Download Reactome Pathway file
#' @description
#'   Download Reactome Pathway file from Reactome web site
#' @param directory Directory/Folder where the data was downloaded.
#' @param mingenenum only pathways consist with larger than 'mingenenum' genes are considered.
#' @param maxgenenum only pathways consist with lower than 'maxgenenum' genes are considered.
#' @param cleanDownload Should downloaded data be deleted?
#' @importFrom downloader download
#' @importFrom AnnotationDbi select
#' @importFrom xfun is_windows
#' @importFrom utils read.delim
#' @import org.Hs.eg.db
#'
#' @return A list of Reactome pathways information
#' @export getReactomePathway
#'
#' @examples
#' \dontrun{
#' reactome <- getReactomePathway(directory = './',
#'                                cleanDownload = T)
#' }
#' @author Yan Li
#'
getReactomePathway <- function(directory,
                               mingenenum = 10,
                               maxgenenum = 200,
                               cleanDownload = FALSE){
  DataDirectory <- paste0(directory, 'NCBI2Reactome.txt')
  fpath <- "https://reactome.org/download/current/NCBI2Reactome.txt"
  if(is_windows()) mode <- "wb" else  mode <- "w"
  message(rep("-",100))
  options(timeout = 10000) # set 10000 second to download the file, default is 60 seconds
  message("o Starting to download Reactome Pathway file")
  if(!file.exists(gsub("\\.txt", "", DataDirectory))){
    download(fpath, DataDirectory, mode = mode)
  }
  message("o Reading Reactome Pathways")
  reactome <- read.delim(DataDirectory, header = F)
  reactome <- reactome[reactome$V6 == 'Homo sapiens', 1:4]
  reactome <- unique(reactome)
  reactome$V1 <- AnnotationDbi::select(org.Hs.eg.db, keys = reactome$V1, columns = 'SYMBOL')$SYMBOL
  reactome <- unique(reactome)
  reactome <- reactome[!(is.na(reactome$V1)),]
  reactome <- reactome[reactome$V4 %in% names(table(reactome$V4)[table(reactome$V4) > mingenenum &
                                                                   table(reactome$V4) < maxgenenum]),]
  if (cleanDownload) unlink(DataDirectory, recursive = T)
  return(reactome)
}
