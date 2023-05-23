#' Generates a dotplot of GSEA results joining all conditions when there are less than 100 results, and an equivalent heatmap in Excel format for all cases 
#' 
#' @param gsea: GSEA results object
#' @param contrast: List of vectors with each contrast to use
#' @param p.adj: Adjusted p-value threshold. Default=0.05
#' @param collectionName: Name of the collection used for the GSEA, e.g. "H" for Hallmark. Default: ""
#' @param resultsDir: Output directory where results will be stored. Default: current working directory

#' 
#' @import ggplot2
#' @import dplyr
#' @import openxlsx
#' 
#' @return Saves a dotplot of GSEA significant results in all conditions, as well as an excel with the same information


makeJoinedDotplot <- function(gsea,contrast,p.adj=0.05,collectionName="", resultsDir=getwd()) {
  
  result_list <- list()
  merged_df <- data.frame()
  # Loop to iterate each GSEA object
  for (i in 1:length(contrast)) {
    
    result <- bind_rows(gsea[[i]]@result[gsea[[i]]@result$p.adjust < p.adj,], .id = "analysis")
    
    if(nrow(results)!=0) { # if there are significant results
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
