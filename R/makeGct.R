#' Create gct file
#' 
#' Creates a gct file from a matrix, to be used in GSEA
#'
#' @param exprMat expression matrix, where rows are identifiers and columns samples
#' @param resultsDir Output directory. Default = working directory
#' @param fileName name of the output file, without extension.
#'
#' @return gct file with the specific structure in specified path
#' @export gctMake

gctMake <-function(exprMat, resultsDir=getwd(), fileName){
  
  rows <- nrow(exprMat)
  cols <- ncol(exprMat)
  
  gct <- file.path(resultsDir, paste(fileName,"gct", sep = "."))
  
  #Generate the headers
  cat("#1.2", '\n', sep="", file = gct)
  cat(rows, '\t', sep="", file = gct, append=TRUE)
  cat(cols, '\n', sep="", file = gct, append=TRUE)
  
  #Generate the columns of data
  colnm <- c("NAME", "Description", colnames(exprMat))
  col1 <- rownames(exprMat)
  col2 <- rep("NA", rows)
  gctfile <- rbind(colnm, data.frame(col1, col2, exprMat, stringsAsFactors = FALSE))
  #Append the data to the file
  write.table(gctfile, file = gct, append = T, sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}
