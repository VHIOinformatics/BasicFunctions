#' Creates an excel file out of a resultsRNAseq file 
#'
#' @param pathRDS: path to RDS file with a data.frame obtained from resultsRNAseq() object
#' @param fileRDS: RDS file with a data.frame obtained from resultsRNAseq() object, without extension
#' @param contrast: List of vectors with each contrast to use
#' @param resultsDir: Output directory
#' @param fileName: Name of the output file, without extension
#' @param pvalue: p-value threshold. Default = NULL, as it is assumed to use an adjusted p-value by default
#' @param padj: adjusted p-value threshold. Default = 0.05
#' @param logFC: abs(logFC) threshold. Default = 1
#' @param add.colors: colors to add if there are more than 20 contrasts. Default = NULL
#'
#' @import openxlsx
#' @import wrapr

#' @return Excel file in resultsDir with fileName
#' @export makeExcelResults

makeExcelResults <- function(pathRDS, fileRDS, contrast,
                                  resultsDir, fileName, pvalue = NULL, 
                                  padj = 0.05, logFC = 1, add.colors = NULL)
  {
  
  message("Reading RDS file ...")
  res<-readRDS(file=file.path(pathRDS, paste0(fileRDS, ".rds")))
  # Change GO columns to eliminate the space at the beginning
  if ("GO.BP" %in% colnames(res)){
    res$GO.BP <- gsub(" ","", res$GO.BP)
  }
  if ("GO.CC" %in% colnames(res)){
    res$GO.CC <- gsub(" ","", res$GO.CC)
  }
  if("GO.MF" %in% colnames(res)){
    res$GO.MF <- gsub(" ","", res$GO.MF)
  }
  #Select heatmap colors columns
  color.values <- res[,grep(colnames(res),pattern=".scaled",fixed = TRUE)]
  
  ## Create a new workbook
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "AllData")
  
  writeData(wb, 1, res)
  
  
  ########### FILTER ############
  
  message("Filter ...")
  for (i in 1:length(contrast)){  # only one comparison
    
    if (is.null(pvalue)) { # pvalue = NULL
      thres.p.type = "adj.P.Val"
      res.filter.p<-res[,grep(colnames(res),pattern=thres.p.type,fixed = TRUE)]
      
      if (!any(res.filter.p < padj)){  # no gene passes the padj threshold, filters by pvalue=0.05
        warning("There are no results lower than the padj threshold. PValue = 0.05 will be used to filter data.") ## IRENE NOTE: may consider asking if they want to increase the padj threshold to 0.1 first?
        p.threshold = 0.05
        thres.p.type = "P.Value"
        res.filter.p<-res[,grep(colnames(res),pattern=thres.p.type,fixed = TRUE)]
      } else { # if there are results passing the padj threshold (res.filter.p < padj)
        p.threshold = padj
        warning(paste("Padj =",p.threshold,"was used to filter data."))
      } 
      
    } else { # pvalue != NULL. pvalue will be provided
      p.threshold = pvalue
      thres.p.type = "P.Value"
      res.filter.p<-res[,grep(colnames(res),pattern=thres.p.type,fixed = TRUE)]  
    }
    
    logFC.col <- res[,grep(colnames(res),pattern="logFC",fixed = TRUE)]
    res.final <- res[which(res.filter.p <= p.threshold & abs(logFC.col) > logFC ),]
    #order by FC columns
    res.f <- res.final[orderv(res.final[paste("FC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")]),]
    # add new sheet
    addWorksheet(wb, sheetName = paste(contrast[[i]][1],"vs",contrast[[i]][2],sep="."))
    writeData(wb, paste(contrast[[i]][1],"vs",contrast[[i]][2],sep=".") , res.f) 
    
  }
  
  #### STYLE ####
  message("Style ...")
  
  sheet.num <- length(contrast)+1
  # For all sheets 
  for (i in 1:sheet.num){
    # Create several styles for columns and rows
    headerStyle1 <- createStyle(fontSize = 10,textDecoration = "Bold",
                                wrapText = TRUE,textRotation = 90,halign = "center")
    
    addStyle(wb, sheet = i, headerStyle1, rows = 1, cols = 1:length(colnames(res)),
             gridExpand = TRUE)
    
    bodyStyle1 <- createStyle(fontSize = 10,
                              wrapText = TRUE, valign = "top", halign = "left")
    addStyle(wb, sheet = i, bodyStyle1, rows = 2:(length(rownames(res))+1),
             cols = 1:length(colnames(res)), gridExpand = TRUE)
    
    headerStyle2 <- createStyle(fontSize = 10, halign = "center",textDecoration = "Bold",
                                wrapText = TRUE, textRotation = 90)
    addStyle(wb, sheet = i, headerStyle2, rows = 1,
             cols = 1:dim(color.values)[2], gridExpand = TRUE)
    
    NumberStyle <- createStyle( fontSize = 10, numFmt = "0.00")
    
    NumberStyle_adj.pval <- createStyle (fontSize = 10, numFmt = "SCIENTIFIC")
    
    HeatmapStyle <- createStyle(fontSize = 10, numFmt = "0")
    addStyle(wb,sheet=i, HeatmapStyle, rows = 2:(nrow(res)+1) , cols = 1:(dim(color.values)[2]),gridExpand=T)
    number.col <- ncol(res) - dim(color.values)[2]
    FCcols <- grep("FC", colnames(res))
    meanCols <- grep("mean", colnames(res))
    adjPvalCols <- grep("adj.P.Val", colnames(res))
    scaleCols <- grep(".scaled", colnames(res))
    
    addStyle(wb, sheet = i, NumberStyle, rows = 2:(nrow(res)+1),
             cols = number.col : ncol(res),gridExpand=T)
    addStyle(wb, sheet = i, NumberStyle, rows = 2:(nrow(res)+1),
             cols = FCcols,gridExpand=T)
    addStyle(wb, sheet = i, NumberStyle, rows = 2:(nrow(res)+1),
             cols = meanCols,gridExpand=T)
    addStyle(wb, sheet = i , NumberStyle_adj.pval, rows = 2:(nrow(res)+1), cols = adjPvalCols,gridExpand=T)
    # Set Heights and Widths
    setRowHeights(wb, sheet = i, rows = 1, heights = 150)
    setRowHeights(wb, sheet = i, rows = 2:(nrow(res)+1), heights = 14)
    
    setColWidths(wb, sheet = i, cols = (ncol(color.values)-1):ncol(res), widths = 8)
    setColWidths(wb, sheet = i, cols = 1:ncol(color.values), widths = 2)
    setColWidths(wb, sheet = i, cols =  number.col : ncol(res), widths = 5)
    setColWidths(wb, sheet = i, cols =  FCcols, widths = 5)
    setColWidths(wb, sheet = i, cols =  meanCols[1]:ncol(res), widths = 4)
    # Fix first row
    freezePane(wb, sheet = i , firstRow = TRUE)
    
    # Change some widths according to specific columns
    if ("AffyID" %in% colnames(res)){
      number.col<-which(colnames(res) == "AffyID")
      setColWidths(wb, sheet = i, cols = number.col, widths = 18)
    }
    if ("Symbol" %in% colnames(res)){
      number.col<-which(colnames(res) == "Symbol")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("Geneid" %in% colnames(res)){
      number.col<-which(colnames(res) == "Geneid")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("mrna" %in% colnames(res)){
      number.col<-which(colnames(res) == "mrna")
      setColWidths(wb, sheet = i, cols = number.col, widths = 4)
    }
    if ("UCSC_symbols" %in% colnames(res)){
      number.col<-which(colnames(res) == "UCSC_symbols")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    
    if ("GO.BP" %in% colnames(res)){
      number.col<-which(colnames(res) == "GO.BP")
      setColWidths(wb, sheet = i, cols = number.col, widths = 12)
    }
    
    if ("GO.CC" %in% colnames(res)){
      number.col<-which(colnames(res) == "GO.CC")
      setColWidths(wb, sheet = i, cols = number.col, widths = 12)
    }
    
    if ("GO.MF" %in% colnames(res)){
      number.col<-which(colnames(res) == "GO.MF")
      setColWidths(wb, sheet = i, cols = number.col, widths = 12)
    }
    
    if ("^path" %in% colnames(res)){
      number.col<-which(colnames(res) == "^path")
      setColWidths(wb, sheet = i, cols = number.col, widths = 4)
    }
    if ("Description" %in% colnames(res)){
      number.col<-which(colnames(res) == "Description")
      setColWidths(wb, sheet = i, cols = number.col, widths = 40)
    }
    if ("Length" %in% colnames(res)){
      number.col<-which(colnames(res) == "Length")
      setColWidths(wb, sheet = i, cols = number.col, widths = 6)
    }
    if ("Strand" %in% colnames(res)){
      number.col<-which(colnames(res) == "Strand")
      setColWidths(wb, sheet = i, cols = number.col, widths = 3)
    }
    if ("Chr" %in% colnames(res)){
      number.col<-which(colnames(res) == "Chr")
      setColWidths(wb, sheet = i, cols = number.col, widths = 5)
    }
    if ("Chrom" %in% colnames(res)){
      number.col<-which(colnames(res) == "Chrom")
      setColWidths(wb, sheet = i, cols = number.col, widths = 4)
    }
    if ("Start" %in% colnames(res)){
      number.col<-which(colnames(res) == "Start")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("Stop" %in% colnames(res)){
      number.col<-which(colnames(res) == "Stop")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    # Heatmap :
    conditionalFormatting( wb,
                           sheet = i,
                           cols = 1:as.numeric(dim(color.values)[2]),
                           rows = 1:as.numeric(dim(color.values)[1] + 1),
                           rule = as.numeric(c(min(color.values),0,max(color.values))),
                           style = c("blue","white", "red"),
                           type = "colorScale"
    )
  }
  
  # COLOURING CONTRASTS: 
  stats<-list()
  #Vector of contrast columns colors
  colors4stats <- c(c("#FFEA00", "#FFC000", "#00B0F0", "#92D050", "#FF6600", "#CCFF99","#CC99FF", "#FF5252", "#5C45FF", "#45FFC7","#fc79f4","#00B0F0", "#9458d1","#c2a03a", "#d1589b","#b3a7cc","#ccf1ff","#1fad66", "#ffeacc", "#f0a1a1" ),add.colors)
  #Only for sheets in contrasts
  for (i in 1:length(contrast)){
    adjp=paste("adj.P.Val",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    p=paste("P.Value",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    fc=paste("FC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    logfc=paste("logFC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    # Select contrast columns
    stats[[i]] <- which(colnames(res) %in% c(adjp,p,fc,logfc))
    # Colouring rows
    conditionalFormatting( wb,
                           sheet = i+1,
                           cols = stats[[i]],
                           rows = 1:(nrow(res)+1),
                           rule = ".",
                           style = createStyle(bgFill = colors4stats[i], fontSize = 10),
                           type = "contains"
    )# Colouring cols
    conditionalFormatting( wb,
                           sheet = i+1,
                           cols = stats[[i]],
                           rows = 1,
                           rule = ".",
                           style = createStyle(bgFill = colors4stats[i], fontSize = 10),
                           type = "contains"
    )
    
  }
  # Only for sheet ALL DATA:
  for (i in 1:length(contrast)){
    adjp=paste("adj.P.Val",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    p=paste("P.Value",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    fc=paste("FC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    logfc=paste("logFC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    # Select contrast columns
    stats[[i]] <- which(colnames(res) %in% c(adjp,p,fc,logfc))
    # Colouring rows
    conditionalFormatting( wb,
                           sheet = 1,
                           cols = stats[[i]],
                           rows = 1:(nrow(res)+1),
                           rule = ".",
                           style = createStyle(bgFill = colors4stats[i], fontSize = 10),
                           type = "contains"
    )# Colouring col
    conditionalFormatting( wb,
                           sheet = 1,
                           cols = stats[[i]],
                           rows = 1,
                           rule = ".",
                           style = createStyle(bgFill = colors4stats[i], fontSize = 10),
                           type = "contains"
    )
    
  }
  
  #### Legend: ####
  message("Legend ...")
  addWorksheet(wb, sheetName = "Legend")
  a<-min(color.values)/5
  b<-max(color.values)/5
  vector.min <- c(min(color.values),a*4,a*3,a*2,a)
  vector.max <- c(b,b*2,b*3,b*4, max(color.values))
  legend.df<-as.numeric(c(vector.min, 0 , vector.max))
  options("openxlsx.borderColour" = "black")
  options("openxlsx.borderStyle" = "thin")
  writeData(wb,  "Legend" , legend.df, borderStyle = getOption("openxlsx.borderStyle", "thin")
            ,startCol = 1,startRow = 1, borders = "columns")
  
  conditionalFormatting( wb,
                         sheet = "Legend",
                         cols = 1,
                         rows = 1:12,
                         rule = as.numeric(c(min(color.values),0,max(color.values))),
                         style = c("blue","white", "red"),
                         type = "colorScale"
  )
  saveWorkbook(wb, file.path(resultsDir,paste(fileName, "xlsx", sep=".")),overwrite = TRUE)
  message("The function was performed successfully")
}



# 09/02/2021: change pvalue of contrast = 1
# change letter size
# 10/02/2021: change widths of columns 

