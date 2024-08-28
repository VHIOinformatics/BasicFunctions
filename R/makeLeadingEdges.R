##' Runs Leading Edges analysis from Gene set enrichment analysis (GSEA)
##' 
##' @param gsea list of GSEA objects (of length n contrasts)
##' @param contrast List with one vector of length 2 for each contrast.
##' @param p.adj Adjusted p-value threshold of significance. Default = 0.05
##' @param collectionName Vector with the name of the collection. It will be appended to the output directory/files name. Example: "H". Default = ""
##' @param gseaDir Directory to save the results
##' @param ... Additional parameters
##' 
##' @export 
##' @import clusterProfiler
##' @import openxlsx 
##' @import dplyr
##' @import tidyr
##' @returns Excel with the genesets in which a significant gene is part of


makeLeadingEdges <- function(gsea, contrast, p.adj=0.05, collectionName="", gseaDir, ...){

  contrast.t <- unlist(lapply(contrast, function(x) paste(x, collapse = ".vs.")))
  contrast.n <- length(contrast)
  #style for all sheets
  headerStyle1 <- createStyle(halign = "left",valign = "center",textDecoration = "Bold", wrapText = TRUE)

  for (i in 1:contrast.n){ # assuming that it was stored in the same order gsea object does not have names for each analysis
    gseaRes <- gsea[[i]]@result
    gseaRes <- gseaRes[gseaRes$p.adjust<p.adj,] # only significant at an adj.p
    a <- strsplit(gseaRes$core_enrichment, split ="/") # list
    names(a) <- gseaRes$ID
    
    # convert the list to dataframe
    dfGenesAll <- stack(a) %>%
      rename(Gene = values, Geneset = ind)
    
    # group the genesets for each gene
    genesGeneSetSum <- dfGenesAll %>%
      group_by(Gene) %>%
      summarise(NGenesets=n(), Genesets = paste(Geneset, collapse = "/")) %>% 
      arrange(desc(NGenesets))
    
    # write to excel
    wb <- createWorkbook()
    addWorksheet(wb, sheetName = "Gene summaries")
    writeData(wb,"Gene summaries", genesGeneSetSum) # POSITIVE NES
    setColWidths(wb, sheet = "Gene summaries", cols = 1:3, widths =c(15,10,100))
    addStyle(wb, sheet = "Gene summaries", headerStyle1, rows = 1, cols = 1:3, gridExpand = TRUE)
    saveWorkbook(wb, file.path(gseaDir, paste("GSEA",collectionName,contrast.t[i],"SignificantGeneSets","GeneSummaries","xlsx",sep=".")), overwrite = TRUE)    

    }


}
