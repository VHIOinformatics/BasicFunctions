#' Create rnk file
#' Creates rnk file using p-values and logFC for sign, to be used in GSEA
#'
#' @param res data frame with columns: Symbol, P.Value and logFC
#' @param colNames Vector with the name of the columns with "Symbol", "P.Value" and "logFC"
#' @param resultsDir Output directory. Default = working directory
#' @param fileName name of the output file, without extension.
#'
#' @return ranked file with the specific structure in specified path
#' @export makeRnk

makeRnk <- function(res, colNames, resultsDir=getwd(), fileName) {
 #add param with contrast to obtain directly from res object
  
  rnk <- file.path(resultsDir, paste(fileName,"rnk", sep = "."))
  
  distrRank <- (-log10(res[, colNames[2]]))*(res[, colNames[3]]/ abs(res[,colNames[3]]))
  resr <- data.frame(Symbol = res[, colNames[1]], distrRank)
  dfbase2 <- aggregate(. ~ Symbol, data = resr, mean)  
  
  write.table(dfbase2, file = rnk, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  
}