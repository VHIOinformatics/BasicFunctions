#' Converts counts matrix to excel format
#'  
#' Creates an excel file from a raw count matrix
#' 
#' @param counts Count matrix.
#' @param resultsDir Directory to save the resulting excel file. Default = working directory
#' @param fileName Name of the excel, without the extension. Default = "rawMatrix"
#' 
#' @import openxlsx
#' @export rawMatrix2Excel
#' 
#' @returns Saves excel file of the count matrix provided.



rawMatrix2Excel <- function(counts, resultsDir=getwd(), fileName="rawMatrix") {
  
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "RawCounts")
  
  writeData(wb, 1, counts, rowNames = TRUE)
  
  headerStyle1 <- createStyle(fontSize = 10,textDecoration = "Bold",
                              wrapText = TRUE,textRotation = 90,halign = "center")
  addStyle(wb, sheet = 1, headerStyle1, rows = 1, cols = 1:(length(colnames(counts))+1),
           gridExpand = TRUE)
  
  bodyStyle1 <- createStyle(fontSize = 10,
                            wrapText = TRUE, valign = "top", halign = "right")
  addStyle(wb, sheet = 1, bodyStyle1, rows = 2:(length(rownames(counts))+1),
           cols = 2:(length(colnames(counts))+1), gridExpand = TRUE)
  
  colStyle1 <- createStyle(fontSize=10,textDecoration="Bold",wrapText=TRUE)
  
  addStyle(wb, sheet = 1, colStyle1, rows = 2:(length(rownames(counts))+1),
           cols = 1, gridExpand = TRUE)
  
  setRowHeights(wb, sheet = 1, rows = 1, heights = 80)
  setColWidths(wb, sheet = 1, cols = 2:(ncol(counts)+1), widths = 6)
  setColWidths(wb, sheet=1, cols=1, widths=10)
  
  freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
  
  saveWorkbook(wb, file.path(resultsDir,paste(fileName, "xlsx", sep=".")),overwrite = TRUE)
}