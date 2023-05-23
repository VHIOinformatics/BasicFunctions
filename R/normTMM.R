#' Function to normalize to TMMs a counts table using edgeR package
#' To be used for visualization, not DEA
#' 
#' @param countsTable table of counts, with genes as rows and samples as columns
#' @param log to be used in function cpm. Default = FALSE
#' 
#' @return Normalized table of counts
#' @import edgeR
#' @export


normTMM <- function(countsTable, log = FALSE, ...){
  
  d <- DGEList(counts=countsTable)
  Norm.Factor <- calcNormFactors(d, method="TMM")
  countsTMM <- cpm(Norm.Factor, log = log, ...)
  
  return(countsTMM)
  
}
  
