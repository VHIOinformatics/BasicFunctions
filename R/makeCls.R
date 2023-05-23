#' Creates a cls file using a vector of conditions, to be used in GSEA
#'
#' @param cond: column of the conditions to be compared of the targets file
#' @param resultsDir: output directory. Default = working directory
#' @param fileName: name of the output file, without extension
#'
#' @return cls file with the specific structure in specified path
#' @export makeCls

makeCls <- function(cond, resultsDir=getwd(), fileName) {
  
  cond <- as.character(cond)
  levels <- unique(cond)
  firstr <- c(length(cond), length(levels),1)
  secndr <- c("#",levels)
  thirdr <- cond
  cls <- file.path(resultsDir, paste(fileName, "cls", sep = "."))
  
  #Create file
  cat(firstr, '\n', file = cls)
  cat(secndr, '\n', file = cls, append = TRUE)
  cat(thirdr, '\n', file = cls, append = TRUE)
  
}
