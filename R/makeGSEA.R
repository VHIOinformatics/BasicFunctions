##' Runs Gene set enrichment analysis (GSEA) using package cluterProfiler 
##'
##' Function that performs a GSEA on a resultsRNAseq object (and a list of contrasts) or a genelist
##' from clusterProfiler package. This function performs an hypergeometric test for each gene set.

##' @param results resultsRNAseq object or a genelist (with names). resultsRNAseq is a data.frame obtained from makeRNAseqResults(). 
##' Must contain a Geneid column and ##' P.Value, adj.P.Val and logFC columns for each contrast (eg. "P.Value.LES_ACT.vs.LES_INACT"). 
##' @param contrast List with one vector of length 2 for each contrast
##' @param gmt Gene sets data.frame as obtained with function clusterProfiler::read.gmt. Must have a column "term" and a column "gene".
##' @param collectionName Vector with the name of the collection. It will be 
##' appended to the output directory/files name (eg. "c5.go.bp")
##' @param resultsDir Character vector with output results directory. Default = working directory.
##' @param minGSSize minimal size of each geneSet for analyzing. Default = 15
##' @param maxGSSize maximal size of genes annotated for testing. Default = 500
##' @param p.value  pvalue cutoff for resultsRNASeq objects. Default = 1, as we want to keep all results
##' @param plots whether to generate plots (bar plot, dot plot, enrichment map and gene concept networks) in the same directory. Default = TRUE
##' @param p.adj threshold of the adjusted p-value to be selected for plots. Default = 0.05
##' @param plotTop Number of maximal gene sets to be plotted in barplot and dotplot. Default = 50

##' 
##' @return Returns a list with enrichment results for each contrast. The list is saved as a RData object 
##' in the resultsDir directory, along with a folder for each contrast containing an 
##' excel file with enrichment results and plots if requested.
##' @export
##' @import clusterProfiler
##' @import openxlsx 
##' @import enrichplot
##' @import gridExtra
##' @import png
##' @import ggplot2
##' @import dplyr
##' @import makeJoinedGSEAplot

makeGSEA <- function(results, contrast, gmt, resultsDir=getwd(), collectionName="", minGSSize=15, maxGSSize=500, p.value=1, plots=TRUE, p.adj=0.05, plotTop=50) {
  
  require(clusterProfiler)
  require(openxlsx)
  require(enrichplot)
  require(gridExtra)
  require(png)
  require(ggplot2)
  require(dplyr)
  
  headerStyle1 <- createStyle(halign = "center",valign = "center",textDecoration = "Bold",
                              wrapText = TRUE) 
  
  gsea=list()
  message("Getting GSEA results...")
  for (i in 1:length(contrast)) {
    # Create a folder for each contrast
    gseaResDir = file.path(resultsDir,paste("GSEA",collectionName,contrast[[i]][1],"vs",contrast[[i]][2], sep="."))
    dir.create(gseaResDir,showWarnings = F)
    # Prepare geneList from results
    p=paste("P.Value",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    logfc=paste("logFC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    
    sign=sign(results[,logfc])
    logP=-log10(results[,p])
    metric=logP*sign
    geneList = metric
    names(geneList) = results$Geneid
    geneList = sort(geneList,decreasing=TRUE)
    
    # GSEA
    set.seed(123)
    gsea[[i]] <- clusterProfiler::GSEA(geneList, TERM2GENE = gmt, verbose=FALSE,seed=T, minGSSize=minGSSize,maxGSSize = maxGSSize, pvalueCutoff = p.value)
    
    
    # save excel with results
    wb <- createWorkbook()
    for (j in 1:length(contrast[[i]])){ # Loop through positive and negative
      # Positive NES
      addWorksheet(wb, sheetName = paste("Enriched in",contrast[[i]][j]))
      if(j==1){
        writeData(wb, j, gsea[[i]]@result[gsea[[i]]@result$NES>0,-c(2,8)])
        } # POSITIVE NES
      if(j==2){
        writeData(wb, j, gsea[[i]]@result[gsea[[i]]@result$NES<0,-c(2,8)])
        } # NEGATIVE NES
      setColWidths(wb, sheet = j, cols = 1:9, widths =c(35,6,10,10,10,10,5,15,25))
      addStyle(wb, sheet = j, headerStyle1, rows = 1, cols = 1:9, gridExpand = TRUE)
      setRowHeights(wb, sheet = j, rows = 1, heights =30)
      
    }
    
    saveWorkbook(wb, file.path(gseaResDir, paste("GSEA",collectionName,contrast[[i]][1],"vs",contrast[[i]][2],"xlsx",sep=".")), overwrite = TRUE)
    
    message(paste0("Results for contrast ",contrast[[i]][1],".vs.",contrast[[i]][2])," saved succesfully")

  }
  save(gsea,file=file.path(resultsDir, paste0("GSEA",collectionName,".RData")))
  
  # make GSEA plots, calling the function
  if (plots) { 
    message("Drawing GSEA plots...")
    makePlotsGSEA(gsea, contrast = contrast, collectionName = collectionName, resultsDir=resultsDir, plotTop=plotTop, p.adj=p.adj)
  }
  
  return(gsea)
  
}



makePlotsGSEA <- function(gsea, contrast, collectionName="", resultsDir=getwd(), plotTop=50, p.adj=0.05){
  
  for (i in 1:length(contrast)) {
    message("Plotting ", paste(contrast[[i]][1],"vs",contrast[[i]][2],sep="."),"\n")
    gseaResDir = file.path(resultsDir,paste("GSEA",collectionName,contrast[[i]][1],"vs",contrast[[i]][2], sep="."))
    dir.create(gseaResDir,showWarnings = F)
    # Create a list with significant [[1]] positive and [[2]] negative results
    gsea.L=list()
    gsea.L[[contrast[[i]][1]]]=gsea[[i]]
    gsea.L[[contrast[[i]][1]]]@result=gsea.L[[contrast[[i]][1]]]@result[gsea.L[[contrast[[i]][1]]]@result$NES>0&gsea.L[[contrast[[i]][1]]]@result$p.adjust<p.adj,]
    gsea.L[[contrast[[i]][2]]]=gsea[[i]]
    gsea.L[[contrast[[i]][2]]]@result=gsea.L[[contrast[[i]][2]]]@result[gsea.L[[contrast[[i]][2]]]@result$NES<0&gsea.L[[contrast[[i]][2]]]@result$p.adjust<p.adj,]
    
    # Generate a general Barplot with NES and pval (maximum to top 30 up and down)
    up=gsea.L[[contrast[[i]][1]]]@result[,c("ID","NES","p.adjust")]
    up=up[1:ifelse(nrow(up)>30,30,nrow(up)),]
    down=gsea.L[[contrast[[i]][2]]]@result[,c("ID","NES","p.adjust")]
    down=down[1:ifelse(nrow(down)>30,30,nrow(down)),]
    data=na.omit(rbind(up,down))
    data$ID=strtrim(data$ID, 70) # maximum label length
    data$ID=factor(data$ID,data$ID[order(data$NES)],ordered=T)
    if(nrow(data)!=0){
      p=ggplot(data, aes(x=ID, y=NES,fill=p.adjust)) + 
        geom_bar(stat = "identity") +
        scale_fill_continuous(low='red', high='blue')+
        ylim(ifelse(min(data$NES)< (-2),min(data$NES),-2),
             ifelse(max(data$NES)> (2),max(data$NES),2)) +
        coord_flip() +
        theme_bw() +
        ggtitle(paste0("GSEA of ", collectionName)) +
        theme(plot.title = element_text(hjust = 0.5,face="bold"))
      
      ggsave(file.path(gseaResDir, 
                       paste0("GSEA.",collectionName,".BarplotNES.",
                              paste0(contrast[[i]][1],"vs",contrast[[i]][2]),".png")),
             plot=p,width = 10,height =8*(ceiling(nrow(data)/80)))
      
    }else{
      png::writePNG(array(0, dim = c(1,1,4)), 
                    file.path(gseaResDir,
                              paste0("GSEA.",collectionName,".BarplotNES.",
                                     paste0(contrast[[i]][1],"vs",contrast[[i]][2]),".png")))
    }  
    for (j in 1:length(contrast[[i]])){ # loop through positive and negative
      if(nrow(gsea.L[[j]]@result)!=0){ # If there are significant results; do plots
        # GSEA Dotplot
        gsea.L[[j]]@result$Description=strtrim(gsea.L[[j]]@result$Description, 70) # maximum label length
        
        p=clusterProfiler::dotplot(gsea.L[[j]], showCategory=40,font.size=8,title= paste0("Enriched in ",names(gsea.L)[j],"\n p.adjust<", p.adj))
        
        ggsave(file.path(gseaResDir, paste0("GSEA.",collectionName,".Dotplot.", names(gsea.L)[j],".png")),
               plot=p,width =9,height =7)
        
        # GSEA Running score
        p=list()
        ## Create GSEA dir for Running Score results
        GSEA_score_dir=paste0(gseaResDir,"/","RunningScore")
        dir.create(GSEA_score_dir)
        ## Set max num. of plots to all DE or to plotTop
        n_plot=ifelse(nrow(gsea.L[[j]]@result)>plotTop,plotTop,nrow(gsea.L[[j]]@result))
        for (u in 1:n_plot){
          p=enrichplot::gseaplot2(gsea.L[[j]],color="green", base_size = 4.5, 
                                  geneSetID = u,title=gsea.L[[j]]$Description[u])
          print(p)
          ggsave(filename = file.path(GSEA_score_dir, paste0(gsea.L[[j]]$Description[u],".RunningScore.",names(gsea.L)[j],".png")),
                 width = 3,height =2.5)
        }
        
        
        # GSEA Gene-Concept networks (top 5 gene sets by default)
        p=clusterProfiler::cnetplot(gsea.L[[j]], foldChange=gsea.L[[j]]@geneList, cex.params=list(category_node=0.7, category_label=0.7, gene_label=0.5),layout = "kk",showCategory = 5)
        p=p+ scale_color_gradient2(name = "-log(p.val)*signFC", low = "blue", mid = "white", high = "red")
        ggsave(file.path(gseaResDir, paste0("GSEA.",collectionName,".GeneConceptNetworks.", 
                                            names(gsea.L)[j], ".png")), plot=p)
        # GSEA EnrichmentMAP (Jaccard index. Plot top 30 by default)
        if(nrow(gsea.L[[j]]@result)>1){
          pt=enrichplot::pairwise_termsim(gsea.L[[j]], method = "JC", semData = NULL, showCategory = 200)
          p <- clusterProfiler::emapplot(pt,cex.params=list(category_label= 0.5),showCategory = 30)
          ggsave(file.path(gseaResDir, paste0("GSEA.",collectionName,".EnrichmentMAP.", 
                                              names(gsea.L)[j], ".png")),plot=p)
        }else{png::writePNG(array(0, dim = c(1,1,4)), file.path(gseaResDir, paste0("GSEA.",collectionName,".EnrichmentMAP.", 
                                                                                   names(gsea.L)[j], ".png")))}
      }else{
        png::writePNG(array(0, dim = c(1,1,4)), file.path(gseaResDir, paste0("GSEA.",collectionName,".Dotplot.", names(gsea.L)[j],".png")))
        png::writePNG(array(0, dim = c(1,1,4)), file.path(gseaResDir, paste0("GSEA.",collectionName,".GeneConceptNetworks.", 
                                                                             names(gsea.L)[j], ".png")))
        png::writePNG(array(0, dim = c(1,1,4)), file.path(gseaResDir, paste0("GSEA.",collectionName,".EnrichmentMAP.", 
                                                                             names(gsea.L)[j], ".png")))
        
      }
    }
  }
  
  makeJoinedDotplot(gsea,contrast,p.adj,collectionName,resultsDir)
  
  message("GSEA performed succesfully!")
  
}



