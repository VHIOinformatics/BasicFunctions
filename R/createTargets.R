#' Creates an empty targets file with the specified names. If a named list of variables is provided, creates the filled targets file.
#'
#' @param names: usually the sample names to be analyzed
#' @param column_list: named list of the columns to be added to the targets file, defining the sample. Usally treatment, timepoint, condition... Default = NULL
#' @param resultsDir: Output directory. Default = current working directory
#' @param fileName: name of the output file, with extension. Default = "targets.txt"
#'
#' @return empty targets file in resultsDir directory, or full targets file if named list of the columns is provided
#' @export createTargets

createTargets <- function(names, column_list=NULL,fileName="targets.txt", resultsDir=getwd()){
  
  if (is.null(column_list)) {
    df <- data.frame("fileName" = names)
  } else {
    df <- data.frame("fileName" = names, column_list)
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
