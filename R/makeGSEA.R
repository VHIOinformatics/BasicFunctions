##' Runs Gene set enrichment analysis (GSEA) using package clusterProfiler 
##'
##' Function that performs a GSEA on a RNAseq.resAnnot object (and a list of contrasts) or a genelist
##' from clusterProfiler package. This function performs an hypergeometric test for each gene set.

##' @param results ResultsRNAseq object or a genelist (with names). resultsRNAseq is a data.frame obtained from `RNAseq.resAnnot()`. 
##' Must contain a Geneid column and ##' P.Value, adj.P.Val and logFC columns for each contrast. Example: "P.Value.LES_ACT.vs.LES_INACT" 
##' @param contrast List with one vector of length 2 for each contrast.
##' @param gmt Gene sets data.frame as obtained with function `clusterProfiler::read.gmt`. Must have a column "term" and a column "gene".
##' @param collectionName Vector with the name of the collection. It will be 
##' appended to the output directory/files name. Example: "c5.go.bp"
##' @param resultsDir Character vector with output results directory. Default = working directory
##' @param minGSSize Minimal size of each geneSet for analyzing. Default = 15
##' @param maxGSSize Maximal size of genes annotated for testing. Default = 500
##' @param p.valuep-value cCutoff for RNAseq.resAnnot objects. Default = 1, as we want to keep all results
##' @param plots Whether to generate plots (bar plot, dot plot, enrichment map and gene concept networks) in the same directory. Default = TRUE
##' @param p.adj Threshold of the adjusted p-value to be selected for plots. Default = 0.05
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


makeGSEA <- function(results, contrast, gmt, resultsDir=getwd(), collectionName="", minGSSize=15, maxGSSize=500, p.value=1, plots=TRUE, p.adj=0.05, plotTop=50) {
  
  
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
    
    message(paste0("\n --- Results for contrast ",contrast[[i]][1],".vs.",contrast[[i]][2])," saved succesfully ---")

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
    message("Plotting ", paste(contrast[[i]][1],"vs",contrast[[i]][2],sep="."))
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



#' Make GSEA dotplot
#' 
#' Generates a dotplot of GSEA results joining all comparisons when there are less than 100 significant genesets results, and an equivalent heatmap in Excel format for all cases 
#' 
#' @param gsea GSEA results object.
#' @param contrast List of vectors with each contrast to use
#' @param p.adj Adjusted p-value threshold. Default = 0.05
#' @param collectionName Name of the collection used for the GSEA, e.g. "H" for Hallmark. Default = ""
#' @param resultsDir Output directory where results will be stored. Default = current working directory

#' 
#' @import ggplot2
#' @import dplyr
#' @import openxlsx
#' @seealso [makeGSEA()]
#' 
#' @return Saves a dotplot of GSEA significant results in all conditions, as well as an excel with the same information
#' @export makeJoinedDotplot


makeJoinedDotplot <- function(gsea,contrast,p.adj=0.05,collectionName="", resultsDir=getwd()) {
  
  result_list <- list()
  merged_df <- data.frame()
  # Loop to iterate each GSEA object
  for (i in 1:length(contrast)) {
    
    result <- bind_rows(gsea[[i]]@result[gsea[[i]]@result$p.adjust < p.adj,], .id = "analysis")
    
    if(nrow(result)!=0) { # if there are significant results
      # add each result to results list
      result_list[[i]] <- result
      result_list[[i]]$analysis <- paste0(contrast[[i]][1],".vs.",contrast[[i]][2])
      df_heat <- result_list[[i]][,c("NES","p.adjust")]
      colnames(df_heat) <- paste0(colnames(df_heat), ".", result_list[[i]]$analysis[[1]])
      
      # create dataframe that will store all results for the later creation of "heatmap excel"
      if (nrow(merged_df)==0) {
        merged_df <- df_heat
      } else {
        merged_df <- merge(merged_df, df_heat, by="row.names", all=TRUE)
        rownames(merged_df) <- merged_df$Row.names
        merged_df$Row.names <- NULL
      }
    }
    
  }
  
  # Add results into a single dataframe
  final_result <- bind_rows(result_list)
  
  # sort the columns of "heatmap excel" dataframe
  merged_df$Gene.Set.ID <- rownames(merged_df)
  res <- merged_df[, sort(colnames(merged_df))]
  
  
  # Remove the collection part from the GS name
  if(strsplit(final_result$ID, "_")[[1]][1] %in% c("HALLMARK","KEGG","REACTOME")){
    collectionName=strsplit(final_result$ID, "_")[[1]][1]
    final_result$ID <- sub("HALLMARK_|KEGG_|REACTOME_","",final_result$ID)
    
  }
  
  # Limit the number of significative results to make a plot.. and the letter size
  if(length(unique(final_result$ID)) <= 100) { # max 100 sig results to make the dotplot
    if (length(unique(final_result$ID)) > 50 & length(unique(final_result$ID)) < 100) {
      charsize=7
    } else {
      charsize=10
    }
    
    # make the Plot
    fplot <- final_result %>% 
      group_by(ID) %>% 
      ggplot(., aes(x = analysis, y = ID, size = p.adjust, col = NES))+ 
      geom_point()+ 
      scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', )+ 
      scale_size(range = c(1,5), trans = 'reverse')+ scale_y_discrete(limits = rev)+ theme_bw()+ ggtitle(paste("Significant GSEA results in",collectionName))+
      theme(axis.title = element_blank(), 
            axis.text = element_text(size = charsize, face = 'bold'), 
            axis.text.x = element_text(angle = 45, hjust =1), 
            plot.title = element_text(size = 18, face = 'bold'), 
            strip.background = element_rect( color="black", fill="black", size=1.5, linetype="solid"), 
            strip.text = element_text(color = 'white', face = 'bold', size = 12))
    
    
    ggsave(file.path(resultsDir, paste0("GSEA.",collectionName,".Total.DotplotNES.png")),
           plot=fplot,width = 10,height =8)
    
  } else { # more than 100 sig results  
    message("Too many significant results. Cannot perform a plot of all as it would be too confusing to interpret. Please look at the equivalent excel file.")
  }
  
  # Create Heatmap Excel with GSEA results for all conditions (NES + pvalue only)
  color.values <- res[,grep(colnames(res),pattern="NES",fixed = TRUE)]
  
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "GSEAresults")
  
  writeData(wb, 1, res[,-grep(colnames(res),pattern="p.adjust",fixed = TRUE)])
  
  headerStyle1 <- createStyle(fontSize = 10,textDecoration = "Bold",
                              wrapText = TRUE,textRotation = 90,halign = "center")
  
  addStyle(wb, sheet = 1, headerStyle1, rows = 1, cols = 1:length(colnames(res)),
           gridExpand = TRUE)
  
  number.col <- ncol(res) - dim(color.values)[2]
  
  # font size depending on pvalue levels
  pvalueLevel1 <- createStyle(fontSize = 8, numFmt = "0.00")
  pvalueLevel2 <- createStyle(fontSize = 10, numFmt = "0.00")
  pvalueLevel3 <- createStyle(fontSize = 12, numFmt = "0.00")
  pvalueLevel4 <- createStyle(fontSize = 14, numFmt = "0.00")
  
  for(i in 1:contrast.n) {
    threshold_p <- paste0("p.adjust.", contrast.t[i])
    column_idx <- match(threshold_p, names(res))
    column_NES <- paste0("NES.",contrast.t[i])
    column_NES_idx <- match(column_NES, names(res))
    
    if (column_NES %in% names(res)) { # if there was no significant results for the comparison the column will not exist, so do it only if it does
      level1Rows <- which(res[,column_idx] < 0.05 & res[,column_idx] >= 0.005)
      level2Rows <- which(res[,column_idx] < 0.005 & res[,column_idx] >= 0.0005)
      level3Rows <- which(res[,column_idx] < 0.0005 & res[,column_idx] >= 0.00005)
      level4Rows <- which(res[,column_idx] < 0.00005)
      addStyle(wb, sheet=1, style=pvalueLevel1, rows = level1Rows+1, cols = column_NES_idx, gridExpand=TRUE)
      addStyle(wb, sheet=1, pvalueLevel2, rows = level2Rows+1, cols = column_NES_idx, gridExpand=TRUE)
      addStyle(wb, sheet=1, pvalueLevel3, rows = level3Rows+1, cols = column_NES_idx, gridExpand=TRUE)
      addStyle(wb, sheet=1, pvalueLevel4, rows = level4Rows+1, cols = column_NES_idx, gridExpand=TRUE)
      
    }
  }
  
  setRowHeights(wb, sheet = 1, rows = 1, heights = 150)
  setColWidths(wb, sheet=1, cols=1, widths=max(nchar(rownames(res)))+1)
  setColWidths(wb, sheet=1, cols=2:number.col, widths=8)
  
  freezePane(wb, sheet = 1 , firstRow = TRUE, firstCol=TRUE)
  
  # Heatmap :
  conditionalFormatting( wb,
                         sheet = 1,
                         cols = 2:as.numeric(dim(color.values)[2] + 1),
                         rows = 2:as.numeric(dim(color.values)[1] + 1),
                         rule = as.numeric(c(min(color.values, na.rm=TRUE),0,max(color.values, na.rm=TRUE))),
                         style = c("blue","white", "red"),
                         type = "colorScale"
  )
  
  # legend of pvalue levels vs fontsize 
  addWorksheet(wb, sheetName = "Legend")
  
  legend.df<-data.frame(p.value = c(5e-2,5e-3,5e-4,5e-5),
                        fontSize = c("NES","NES","NES","NES"))
  
  options("openxlsx.borderColour" = "black")
  options("openxlsx.borderStyle" = "thin")
  writeData(wb,  "Legend" , legend.df, borderStyle = getOption("openxlsx.borderStyle", "thin")
            ,startCol = 1,startRow = 1, borders = "columns")
  
  addStyle(wb, sheet=2, pvalueLevel1, rows = 2, cols = 2, gridExpand=TRUE)
  addStyle(wb, sheet=2, pvalueLevel2, rows = 3, cols = 2, gridExpand=TRUE)
  addStyle(wb, sheet=2, pvalueLevel3, rows = 4, cols = 2, gridExpand=TRUE)
  addStyle(wb, sheet=2, pvalueLevel4, rows = 5, cols = 2, gridExpand=TRUE)
  
  
  excelFileName="GSEA.Results_HeatmapNES"
  saveWorkbook(wb, file.path(resultsDir,paste(excelFileName, collectionName, "xlsx", sep=".")),overwrite = TRUE)
  
}

