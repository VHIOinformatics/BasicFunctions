##' Runs over representation analysis using enrichr function of cluterProfiler 
##'
##' Function that performs a gene-set enrichment analysis using the function enricher
##' on a resultsRNAseq object (and a list of contrasts) or a genelist
##' from clusterProfiler package. This function performs an hypergeometric test for each gene set.
##' 
##' @param results resultsRNAseq object or a genelist (with names). resultsRNAseq is a data.frame obtained from resultsRNAseq(). 
##' Must contain a Geneid column and ##' P.Value, adj.P.Val and logFC columns for each contrast (eg. "P.Value.LES_ACT.vs.LES_INACT"). 
##' ORA will be performed on the genelist if contrast parameter is NULL
##' @param contrast List with one vector of length 2 for each contrast. Default = NULL
##' @param gmt Gene sets data.frame as obtained with function clusterProfiler::read.gmt. Must have a column "term" and a column "gene".
##' @param collectionName Vector with the name of the collection. It will be 
##' appended to the output directory/files name (eg. "c5.go.bp")
##' @param resultsDir Character vector with output results directory. Default is working directory.
##' @param minGSSize minimal size of each geneSet for analyzing. Default = 15
##' @param maxGSSize maximal size of genes annotated for testing. Default = 500
##' @param p.value  pvalue cutoff for resultsRNASeq objects. Default = 0.05
##' @param p.adjust adjusted pvalue cutoff for resultsRNASeq objects. Default = 0.05
##' @param logFC logFC cutoff for resultsRNASeq objects. Default = 0
##' @param plots whether to generate plots (bar plot, dot plot, enrichment map and gene concept networks) in the same directory. Default = TRUE
##' @param plotTop Number of maximal gene sets to be plotted in barplot and dotplot. Default = 50
##' @param plotP.adjust threshold of the adjusted p-value to be selected for plots. Default = 0.05
##' 
##' @return Returns a list with enrichment results for each contrast. The list is saved as a RData object 
##' in the resultsDir directory, along with a folder for each contrast containing an 
##' excel file with enrichment results.

##' @export
##' @import clusterProfiler
##' @import openxlsx 
##' @import plotORA

makeORA <-function(results, contrast = NULL, gmt, collectionName = "", resultsDir = getwd(), minGSSize = 15, maxGSSize = 500, p.value = 0.05, p.adjust = 0.05, logFC=0, plots = TRUE, plotTop = 50, plotP.adjust = 0.05) {
  
  # If resultsRNAseq and contrasts, will run per contrast and separate UP and DOWN
  require(clusterProfiler)
  require(openxlsx)
  
  headerStyle1 <- createStyle(halign = "center", valign = "center", 
                              textDecoration = "Bold", wrapText = TRUE)

  enrichment = list()
  
  if (is.null(contrast)) {
    
        for (i in 1:length(results)) {
          message("Running enrichment ",i,": ",names(results)[i],"\n")
          
          ORAresultsDir = file.path(resultsDir, paste("Enrichment", collectionName, 
                                            names(results)[i], sep = "."))
          dir.create(ORAresultsDir, showWarnings = F)
          
          
          enrichment[[i]] <- clusterProfiler::enricher(results[[i]], TERM2GENE = gmt,    minGSSize = minGSSize, maxGSSize = maxGSSize, pvalueCutoff = 1)
          
          if (plots) {
            plotORA(enrichment[[i]], plotName = names(results)[i], collectionName = collectionName, resultsDir = ORAresultsDir, plotTop = 50, plotP.adjust = plotP.adjust)
          }
          
          names(enrichment)[i] <- names(results)[i]
          wb <- createWorkbook()
          addWorksheet(wb, sheetName = paste0("Enrichment ", collectionName))
          writeData(wb, 1, enrichment[[i]]@result[-c(2,9)])
          setColWidths(wb, sheet = 1, cols = 1:7, widths = c(35, 6, 10, 10, 10, 10, 25))
          addStyle(wb, sheet = 1, headerStyle1, rows = 1, cols = 1:7, gridExpand = TRUE)
          setRowHeights(wb, sheet = 1, rows = 1, heights = 30)
          saveWorkbook(wb, file.path(ORAresultsDir, paste("Enrichment", collectionName, names(results)[i], "xlsx", sep = ".")), overwrite = TRUE)
        }
    
  } else {
    
        for (i in 1:length(contrast)) {
          enrichment[[i]] = list()
          contrastName <- paste(contrast[[i]][1], "vs", contrast[[i]][2], sep = ".")
          names(enrichment)[i] = contrastName
          message("Running enrichment ",i,": ",contrastName,"\n")
          
          ORAresultsDir = file.path(resultsDir, paste("Enrichment", collectionName, contrastName, sep = "."))
          dir.create(ORAresultsDir, showWarnings = F)
          
          adjp = paste("adj.P.Val", contrastName, sep = ".")
          p = paste("P.Value", contrastName, sep = ".")
          logfc = paste("logFC", contrastName, sep = ".")
          
          #Select lists of genes
          genes_UP <- results[results[,p] < p.value & results[,adjp] < p.adjust & results[,logfc] > logFC,"Geneid" ]
          genes_DOWN <- results[results[,p] < p.value & results[,adjp] < p.adjust & results[,logfc] < (-logFC),"Geneid" ]
          
          geneList = list(genes_UP, genes_DOWN)
          names(geneList) = contrast[[i]]
          
          enrichment[[i]][[names(geneList)[1]]] <- clusterProfiler::enricher(geneList[[1]], TERM2GENE = gmt, universe = results$Geneid, minGSSize = minGSSize, maxGSSize = maxGSSize, pvalueCutoff = 1)
          enrichment[[i]][[names(geneList)[2]]] <- clusterProfiler::enricher(geneList[[2]], TERM2GENE = gmt, universe = results$Geneid, minGSSize = minGSSize, maxGSSize = maxGSSize, pvalueCutoff = 1)
          
          if (plots) {
            
            plotORA(enrichment[[i]][[1]], plotName = names(geneList)[1], collectionName = collectionName, resultsDir = ORAresultsDir, plotTop = 50, plotP.adjust = 0.05)
            plotORA(enrichment[[i]][[2]], plotName = names(geneList)[2], collectionName = collectionName, resultsDir = ORAresultsDir, plotTop = 50, plotP.adjust = 0.05)
            
          }
          
          wb <- createWorkbook()
          for (j in 1:length(contrast[[i]])) {
            
            addWorksheet(wb, sheetName = paste("Enriched in", contrast[[i]][j]))
            writeData(wb, j, enrichment[[i]][[names(geneList)[j]]]@result[-c(2,9)])
            
            setColWidths(wb, sheet = j, cols = 1:7, widths = c(35, 6, 10, 10, 10, 10, 25))
            addStyle(wb, sheet = j, headerStyle1, rows = 1, cols = 1:7, gridExpand = TRUE)
            setRowHeights(wb, sheet = j, rows = 1, heights = 30)
            
          }
          
          saveWorkbook(wb, file.path(ORAresultsDir, paste("Enrichment", collectionName, contrastName, "xlsx", sep = ".")), overwrite = TRUE)
        }
  }
  
  save(enrichment, file = file.path(resultsDir, paste("Enrichment", collectionName, "RData", sep=".")))
  
  message(i," enrichment done!\n")
  return(enrichment)
  
  
}


