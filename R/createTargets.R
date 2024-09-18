#' Create target file
#' 
#' Creates an empty targets file with the specified names. If a named list of variables is provided, creates the filled targets file.
#'
#' @param names Vector of the sample names to be analyzed. Usually the column names of the count table
#' @param column_list Named list of the columns to be added to the targets file, defining the sample. Usually treatment, timepoint, condition... Default = NULL
#' @param resultsDir Output directory. Default = current working directory
#' @param fileName Name of the output file, with the extension. Default = "targets.txt"
#'
#' @return Saves a tabulated table at resultsDir, with the columns of the conditions are provided as a named list, or an empty file if nothing is provided (to be filled manually later).
#' @export createTargets
#' @examples
#' # prepare columns
#' treat <- c("Untreat","Untreat","Untreat","Treat","Treat","Treat") 
#' weeks <- c("4w","9w","13w","4w","9w","13w")
#' cond <- paste(treat,weeks, sep=".")
#' # create list of the columns, must be with names for proper conversion to column name
#' column_list <- list(Treat=treat,Weeks=weeks,Cond=cond)
#' 
#' createTargets(names=colnames(countsM), column_list = column_list, resultsDir=resultsDir)
#' 
#' # after the file has created, we read it as table
#' targets <- read.table(file.path(resultsDir,"targets.txt"), dec=",", sep="\t", header=T)

createTargets <- function(names, column_list=NULL,fileName="targets.txt", resultsDir=getwd()){
  
  if (is.null(column_list)) {
    df <- data.frame("sampleName" = names)
  } else {
    df <- data.frame("sampleName" = names, column_list)
  }
  
  if(file.exists(file.path(resultsDir,fileName))) {
    message("Careful! This targets file already exists.")
    overwrite <- readline(prompt="Do you want to overwrite it? (yes/no) ")
    if (overwrite=="yes") {
      fileName = fileName
    } else if (overwrite=="no") {
      newFilename <- readline(prompt="Please, enter a new name for the targets file (including extension): ")
      fileName = newFilename
    }
    
  }
  
  write.table(df, file = file.path(resultsDir, fileName),
              dec = ",", sep = "\t", row.names = FALSE, quote = FALSE)
  
}
