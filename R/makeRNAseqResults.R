#' Creates results of RNAseq analysis, including GO annotations 
#' 
#' @param exprMat Matrix of expression in CPM(E expression provided by voom), TMM or in counts, with rownames = Geneid
#' @param annotMat Matrix with annotations obtained with `getAnnot` or `makeAnnot` functions
#' @param cond Vector with sample conditions used in the contrasts, same order as colnames(exprMat)
#' @param fitMain Object obtained after contrasts.fit() function from limma package is applied. Default = fit.main
#' @param contrast List of vectors with each contrast to use
#' @param species species used to retrieve annotations (human or mouse). Default = "human"
#' @param GO Annotation with GO. Default = TRUE
#' @param geneidCol Column of annotMat that contains gene identifier
#' @param idType type of gene identifier in column geneidCol. Possible values are "SYMBOL" or "ENSEMBL". Default = "SYMBOL"
#' @param resultsDir Output directory where results will be stored.
#' @param resAnnotFilename Filename of the csv and rds files that will be saved, without the extension. Default = "resultsAnnot"
#' @param Excel Create an excel of the results. Default = TRUE
#' @param excelFilename Name of the excel file, without the extension
#' @param pvalue p-value threshold. Default = NULL, as it is assumed to use an adjusted p-value by default
#' @param padj: adjusted p-value threshold. Default = 0.05
#' @param logFC: abs(logFC) threshold. Default = 1
#' @param add.colors: colors to add if there are more than 20 contrasts. Default = NULL
#'
#' @import gtools
#' @import GO.db
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import makeExcelResults
#' 
#' @returns Results of RNA-seq analysis. Returns results object, and saves outputs in csv and rds format to resultsDir, as well as an excel
#' 
#' @export RNAseq.resAnnot  

# library(gtools)
# library(GO.db)
# library(org.Hs.eg.db)
# library(org.Mm.eg.db)

RNAseq.resAnnot <- function(exprMat, annotMat, cond, fitMain = fit.main, contrast, species="human", GO=TRUE, geneidCol = "Geneid", idType="SYMBOL", resultsDir, resAnnotFilename="resultsAnnot", Excel=TRUE, excelFilename, pvalue = NULL, padj = 0.05, logFC = 1, add.colors = NULL) {
  
  #Obtain contrasts with limma
  ConList <- vector("list", length(contrast)) 
  for (i in 1:length(contrast)) {
    ConList[[i]] <- topTable(fit.main,n=Inf,coef=i, adjust="fdr")[,c("logFC","P.Value","adj.P.Val")]
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
    
    GENENAME <- select(ann.database, keys = annotMat.s[,geneidCol], columns = c("GENENAME"), keytype = idType)
    GO <- select(ann.database, keys = annotMat.s[, geneidCol], columns = c("GO"), keytype = idType)
    PATH <- select(ann.database, keys=annotMat.s[,geneidCol], columns=c("PATH"), keytype = idType)
    
    GENENAME.agg <-aggregate(GENENAME, by=list(GENENAME[,idType]), FUN=function(x) paste(x, collapse="//"))
    GENENAME.agg <- GENENAME.agg[, c("Group.1", "GENENAME")]
    GENENAME.agg.s <- GENENAME.agg[order(GENENAME.agg$Group.1),]
    
    GO.Term <- select(GO.db, keys=GO$GO, columns=c("TERM"), keytype="GOID")
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
    GO.annot.all <- data.frame(SYMBOL=GO.annot[,idType], GO.BP, GO.CC, GO.MF)
    
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
    
    GO.annot.desg <- data.frame(SYMBOL=symb.vect, GO.BP.p, GO.CC.p, GO.MF.p)
    GO.annot.agg.s <- GO.annot.desg[order(as.character(GO.annot.desg[,idType])),]
    
    annotMatNEW <- cbind(annotMat.s[,c(1:ncol(annotMat.s))],
                         GENENAME.agg.s$GENENAME,
                         GO.annot.agg.s[,c(2:ncol(GO.annot.agg.s))])
    
    colnames(annotMatNEW) <- c(colnames(annotMat.s[,c(1:ncol(annotMat.s))]),
                               "Description", "GO.BP", "GO.CC", "GO.MF")
    
  } else {
    #Sort annotation matrix
    annotMatNEW <- annotMat[match(rownames(res_scaled),annotMat$Geneid),]
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
