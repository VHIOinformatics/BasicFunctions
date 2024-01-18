#' Obtain results of RNAseq analysis
#' 
#' Creates results of differential expression RNAseq analysis, which may including GO annotations and an excel with the total and per comparison results 
#' 
#' @param exprMat Matrix of expression in CPM(E expression provided by voom), TMM or in counts, with rownames = Geneid.
#' @param annotMat Matrix with annotations obtained with `getAnnot` or `makeAnnot` functions.
#' @param cond Vector with sample conditions used in the contrasts, same order as colnames(exprMat).
#' @param fitMain Object obtained after contrasts.fit() function from limma package is applied. Default = fit.main
#' @param contrast List of vectors with each contrast to use.
#' @param species Species used to retrieve annotations (human or mouse). Default = "human"
#' @param GO Add annotation with GO. Default = TRUE
#' @param geneidCol Column of annotMat that contains gene identifier. Default = "Geneid"
#' @param idType type of gene identifier in column geneidCol. Possible values are "SYMBOL" or "ENSEMBL". Default = "SYMBOL"
#' @param resultsDir Output directory where results will be stored. Default = current working directory
#' @param resAnnotFilename Filename of the csv and rds files that will be saved, without the extension. Default = "resultsAnnot"
#' @param Excel Create an excel of the results. Default = TRUE
#' @param excelFilename Name of the excel file, without the extension
#' @param pvalue p-value threshold. Default = NULL, as it is assumed to use an adjusted p-value by default
#' @param padj Adjusted p-value threshold. Default = 0.05
#' @param logFC Abs(logFC) threshold. Default = 1
#' @param add.colors Colors to add if there are more than 20 contrasts. Default = NULL
#'
#' @import gtools
#' @import GO.db
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import AnnotationDbi
#' @seealso [makeExcelResults(), getAnnot(), makeAnnot()]
#' 
#' @returns Results of RNA-seq analysis. Returns results object, and saves outputs in csv and rds format to resultsDir. By default, results are also saved in Excel format.
#' 
#' @export RNAseq.resAnnot 


RNAseq.resAnnot <- function(exprMat, annotMat, cond, fitMain = fit.main, contrast, species="human", GO=TRUE, geneidCol = "Geneid", idType="SYMBOL", resultsDir=getwd(), resAnnotFilename="resultsAnnot", Excel=TRUE, excelFilename, pvalue = NULL, padj = 0.05, logFC = 1, add.colors = NULL) {
  
  #Obtain contrasts with limma
  ConList <- vector("list", length(contrast)) 
  for (i in 1:length(contrast)) {
    ConList[[i]] <- topTable(fitMain,n=Inf,coef=i, adjust="fdr")[,c("logFC","P.Value","adj.P.Val")]
    ConList[[i]] <- ConList[[i]][order(rownames(ConList[[i]])),]
  }
  
  #Scale values from count matrix (normalized or not) for the heatmap
  res<-exprMat[order(rownames(exprMat)),]
  res_centered<-res-apply(res,1,mean)
  res_scaled<-res_centered/apply(res,1,sd)
  colnames(res_scaled) <- paste(colnames(res_scaled),"scaled",sep=".")
  
  #Build matrix with mean expression per condition
  u.cond <- unique(cond)
  mean.matrix <- t(apply(res,1, function(x) tapply(x,cond,mean)))
  mean.matrix <- mean.matrix[,u.cond] #proper order
  colnames(mean.matrix) <- paste("mean",u.cond, sep=".")
  
  #Compute FC using previously computed mean values
  FC.matrix <- matrix(data= NA, nrow=nrow(exprMat), ncol=length(contrast))
  col.FC.Names <- vector()
  for (ic in 1:length(contrast)) {
    FC.matrix[,ic] <- 2^abs(ConList[[ic]]$logFC) * sign(ConList[[ic]]$logFC)
    col.FC.Names <- c(col.FC.Names, paste("FC",contrast[[ic]][1], "vs", contrast[[ic]][2], sep="."))
  }
  colnames(FC.matrix) <- col.FC.Names
  
  #Build topDiff matrix with c("FC", "logFC","P.Value","adj.P.Val")
  topDiff.mat <- matrix(data= NA, nrow=nrow(exprMat), ncol=length(contrast)*4)
  col.topDiff.Names <- vector()
  for (ic in 1:length(contrast)) {
    icc <- 1+(4*(ic-1))
    cont.name <- paste(contrast[[ic]][1], "vs", contrast[[ic]][2], sep=".")
    topDiff.mat[,icc] <- FC.matrix[,ic]
    topDiff.mat[,icc+1] <- ConList[[ic]][,1]
    topDiff.mat[,icc+2] <- ConList[[ic]][,2]
    topDiff.mat[,icc+3] <- ConList[[ic]][,3]
    col.topDiff.Names <- c(col.topDiff.Names, col.FC.Names[ic],
                           paste(colnames(ConList[[ic]]), cont.name, sep="."))
  }
  colnames(topDiff.mat) <- col.topDiff.Names
  
  
  #If GO annotations are requested, call function to annotate
  if(GO) {
    message("Adding GO annotations...")
    annotMat.s <- annotMat[order(annotMat[,geneidCol]),]
    #library(KEGG.db)
    #Sort and build matrix with annotations
    if (species =="human") {
      ann.database <- org.Hs.eg.db
    } else if (species=="mouse") {
      ann.database <- org.Mm.eg.db
    }
    
    GENENAME <- AnnotationDbi::select(ann.database, keys = annotMat.s[,geneidCol], columns = c("GENENAME"), keytype = idType)
    GO <- AnnotationDbi::select(ann.database, keys = annotMat.s[, geneidCol], columns = c("GO"), keytype = idType)
    PATH <- AnnotationDbi::select(ann.database, keys=annotMat.s[,geneidCol], columns=c("PATH"), keytype = idType)
    
    GENENAME.agg <-aggregate(GENENAME, by=list(GENENAME[,idType]), FUN=function(x) paste(x, collapse="//"))
    GENENAME.agg <- GENENAME.agg[, c("Group.1", "GENENAME")]
    GENENAME.agg.s <- GENENAME.agg[order(GENENAME.agg$Group.1),]
    
    GO.Term <- AnnotationDbi::select(GO.db, keys=GO$GO, columns=c("TERM"), keytype="GOID")
    GO.annot <- cbind(GO, GO.Term)
    
    GO.BP <- vector(mode="character", length=nrow(GO.annot))
    GO.CC <- vector(mode="character", length=nrow(GO.annot))
    GO.MF <- vector(mode="character", length=nrow(GO.annot))
    
    for (r in c(1:nrow(GO.annot))) {
      Ont <- GO.annot$ONTOLOGY[r]
      if (!is.na(Ont)){
        if(Ont == "BP"){
          GO.BP[r] <- paste(GO.annot[r,c("ONTOLOGY","GOID", "TERM")], collapse=" ")
        } else if (Ont == "CC") {
          GO.CC[r] <- paste(GO.annot[r,c("ONTOLOGY","GOID", "TERM")], collapse=" ")
        } else if (Ont == "MF") {
          GO.MF[r] <- paste(GO.annot[r,c("ONTOLOGY","GOID", "TERM")], collapse=" ")
        }
      }
    }
    if(idType=="SYMBOL") {
      GO.annot.all <- data.frame(SYMBOL=GO.annot[,idType], GO.BP, GO.CC, GO.MF)
    } else if (idType=="ENSEMBL") {
      GO.annot.all <- data.frame(ENSEMBL=GO.annot[,idType], GO.BP, GO.CC, GO.MF)
    }
    
    symb.vect=unique(GO.annot[,idType])
    GO.BP.p <- vector(mode="character", length=length(symb.vect))
    GO.CC.p <- vector(mode="character", length=length(symb.vect))
    GO.MF.p <- vector(mode="character", length=length(symb.vect))
    
    for(si in  c(1:length(symb.vect))) {
      symb <- symb.vect[si]
      GO.mat <- GO.annot.all[GO.annot.all[,idType] == symb,]
      GO.BP.p[si] <- paste(unique(GO.mat$GO.BP), collapse="")
      GO.CC.p[si] <- paste(unique(GO.mat$GO.CC), collapse="")
      GO.MF.p[si] <- paste(unique(GO.mat$GO.MF), collapse="")
      
      GO.BP.p[si] <- gsub("BP", "//", GO.BP.p[si])
      GO.CC.p[si] <- gsub("CC", "//", GO.CC.p[si])
      GO.MF.p[si] <- gsub("MF", "//", GO.MF.p[si])
      
      GO.BP.p[si] <- sub("//", "", GO.BP.p[si])
      GO.CC.p[si] <- sub("//", "", GO.CC.p[si])
      GO.MF.p[si] <- sub("//", "", GO.MF.p[si])
      
    }
    
    if(idType=="SYMBOL") {
      GO.annot.desg <- data.frame(SYMBOL=symb.vect, GO.BP.p, GO.CC.p, GO.MF.p)
    } else if(idType=="ENSEMBL") {
      GO.annot.desg <- data.frame(ENSEMBL=symb.vect, GO.BP.p, GO.CC.p, GO.MF.p)
    }
    GO.annot.agg.s <- GO.annot.desg[order(as.character(GO.annot.desg[,idType])),]
    
    annotMatNEW <- cbind(annotMat.s[,c(1:ncol(annotMat.s))],
                         GENENAME.agg.s$GENENAME,
                         GO.annot.agg.s[,c(2:ncol(GO.annot.agg.s))])
    
    colnames(annotMatNEW) <- c(colnames(annotMat.s[,c(1:ncol(annotMat.s))]),
                               "Description", "GO.BP", "GO.CC", "GO.MF")
    
  } else {
    #Sort annotation matrix
    annotMatNEW <- annotMat[match(rownames(res_scaled),annotMat[,geneidCol]),]
  }
  
  #Make sure everything is in same order before merging and returning
  a=all.equal(rownames(res_scaled), rownames(res))
  b=all.equal(rownames(res_scaled), rownames(ConList[[1]]))
  c=all.equal(rownames(res_scaled), rownames(mean.matrix))
  d=all.equal(rownames(res_scaled), annotMatNEW$Geneid)
  tryCatch(
    expr = {
      a&b&c&d #If all objects are in same order
      message("All objects are in same order. Results dataframe successfully created!")
      
      #Build final matrix only if above expression is true
      res.annot<-data.frame(res_scaled, annotMatNEW,
                            topDiff.mat, mean.matrix,
                            res)
      
      
      #return(res.annot)
      
    },
    error = function(e){ #If  a&b&c&d generated error, pop up following message:
      message('Something is not in the correct order. Please check that all inputs have same length and Geneid:')
      
      message('all.equal(rownames(res_scaled), rownames(res): ',a)
      message('all.equal(rownames(res_scaled), rownames(ConList[[1]])): ',b)
      message('all.equal(rownames(res_scaled), rownames(mean.matrix)): ',c)
      message('all.equal(rownames(res_scaled), annotMatNEW$Geneid): ',d)
      
    }
  ) 
  
  message("Saving results...")
  csv.file = paste0(resAnnotFilename,".csv")
  rds.file = paste0(resAnnotFilename,".rds")
  if(file.exists(file.path(resultsDir,csv.file)) & file.exists(file.path(resultsDir,rds.file))) {
    message(paste0("Csv and rds files already exist with the name '",resAnnotFilename,"' at ",resultsDir))
    overwrite <- readline(prompt="Do you want to overwrite? (yes/no) ")
    if (overwrite=="yes") {
      csv.file = csv.file
      csv.file = csv.file
    } else if (overwrite=="no") {
      newResAnnotFilename <- readline(prompt="Please, enter a new name that doesn't exist: ")
      csv.file = paste0(newResAnnotFilename,".csv")
      rds.file = paste0(newResAnnotFilename,".rds")
      
    } else {
      message("Answer not valid.")
      newResAnnotFilename <- readline(prompt="Please, enter a new name for the results files: ")
      csv.file = paste0(newResAnnotFilename,".csv")
      rds.file = paste0(newResAnnotFilename,".rds")
    }
  }
  
  write.csv2(res.annot,file = file.path(resultsDir,csv.file), row.names = FALSE)
  saveRDS(res.annot, file = file.path(resultsDir, rds.file))
  
  message(paste("Table of results was saved at:",resultsDir))
  
  
  if(Excel) {
    message("Creating excel with results...")
    makeExcelResults(resultsDir, resAnnotFilename, contrast, resultsDir, excelFilename, pvalue, padj, logFC, add.colors)
  } else {
    message("Excel not created. If you want to create it later, you can use 'makeExcelResults()' function.")
  }
  
  return(res.annot)
  
}


#' Create RNAseq results in Excel
#' 
#' Creates an excel file out of a ResultsRNAseq file 
#'
#' @param pathRDS path to RDS file with a data.frame obtained from `RNAseq.resAnnot()` object
#' @param fileRDS RDS file with a data.frame obtained from resultsRNAseq() object, without extension
#' @param contrast List of vectors with each contrast to use
#' @param resultsDir Output directory
#' @param fileName Name of the output file, without extension
#' @param pvalue p-value threshold. Default = NULL, as it is assumed to use an adjusted p-value by default
#' @param padj adjusted p-value threshold. Default = 0.05
#' @param logFC abs(logFC) threshold. Default = 1
#' @param add.colors colors to add if there are more than 20 contrasts. Default = NULL
#'
#' @import openxlsx
#' @import wrapr
#' @seealso [RNAseq.resAnnot()] 
#'
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
    if (length(contrast) == 1) {
      res.final <- res[which(res.filter.p <= p.threshold & abs(logFC.col) > logFC ),]
    } else {
      res.final <- res[which(res.filter.p[i] <= p.threshold & abs(logFC.col[i]) > logFC ),]
    }
    
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
      if (grepl("^ENS",res$Geneid[1])) { #if ensembl id make it bigger
        setColWidths(wb, sheet = i, cols = number.col, widths = 16)
      } else {
        setColWidths(wb, sheet = i, cols = number.col, widths = 8)
      }
      
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
  message("Saving results ...")
  saveWorkbook(wb, file.path(resultsDir,paste(fileName, "xlsx", sep=".")),overwrite = TRUE)
  message("The function was performed successfully")
}



# 09/02/2021: change pvalue of contrast = 1
# change letter size
# 10/02/2021: change widths of columns 
