##' Runs Gene set enrichment analysis (GSEA) using package clusterProfiler 
##'
##' Function that performs a GSEA on a RNAseq.resAnnot object (and a list of contrasts) or a genelist
##' from clusterProfiler package. This function performs an hypergeometric test for each gene set.

##' @param results ResultsRNAseq object, a data.frame obtained from `RNAseq.resAnnot()`. 
##' Must contain a Geneid column and ##' P.Value, adj.P.Val and logFC columns for each contrast. Example: "P.Value.LES_ACT.vs.LES_INACT" 
##' @param contrast List with one vector of length 2 for each contrast.
##' @param gmt Gene sets data.frame as obtained with function `clusterProfiler::read.gmt`. Must have a column "term" and a column "gene".
##' @param collectionName Vector with the name of the collection. It will be 
##' appended to the output directory/files name. Example: "c5.go.bp". Default = ""
##' @param resultsDir Character vector with output results directory. Default = working directory
##' @param minGSSize Minimal size of each geneSet for analyzing. Default = 15
##' @param maxGSSize Maximal size of genes annotated for testing. Default = 500
##' @param p.value p-value cutoff for RNAseq.resAnnot objects. Default = 1, as we want to keep all results
##' @param p.adj Threshold of the adjusted p-value to be selected for plots. Default = 0.05
##' @param rankMethod By default, the ranking is performed with the log p-value and the sign of fold change ("sign"). However, ranking by both log p-value and logFC is possible when specifying "logFC". Default = "sign"
##' @param plots Whether to generate plots (bar plot, dot plot, enrichment map and gene concept networks) in the same directory. Default = TRUE
##' @param plots2Make List of plots that will be performed by makeGSEAplots, if plots = TRUE. By default, all of them will be performed. Default =  c("Barplot", "Dotplot", "RunScoreSig", "GCNetwork", "EnrichMAP")
##' @param plotTop Number of maximum gene sets to be plotted in barplot and dotplot, and enrichment map. Default = 30
##' @param plotTopRS Number of maximum plots to perform the running score. Default = 50
##' @param JoinedDotplot Whether to perform a joined dotplot that compares the overlap of the enriched terms between comparisons. Default = TRUE
##' @param leadingEdges Whether to make an excel including how many and which gene sets a significant gene is part of. Will make a folder "leadingEdges" with the excels per comparison. Default = FALSE
##' @param contrast.order Order in which the contrasts will be shown in the JoinedDotplot. If NULL, it will use the order as shown in contrast. Default = NULL
##' @param ... Additional parameters

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
##' @import stringr


makeGSEA <- function(results, contrast, gmt, resultsDir = getwd(), collectionName = "", minGSSize = 15, maxGSSize = 500, p.value = 1, p.adj = 0.05, rankMethod = "sign", plots = TRUE, plots2Make = c("Barplot", "Dotplot", "RunScoreSig", "GCNetwork", "EnrichMAP"), plotTop = 30, plotTopRS = 50, JoinedDotplot = TRUE, contrast.order = NULL, leadingEdges = FALSE, ...) {
  
  
  gsea = list()
  message("Getting GSEA results...")
  for (i in 1:length(contrast)) {
    # Create a folder for each contrast
    gseaResDir = file.path(resultsDir, paste("GSEA", collectionName, contrast[[i]][1], "vs", contrast[[i]][2], sep = "."))
    dir.create(gseaResDir, showWarnings = F)
    
    # Prepare geneList from results
    pCol = paste("P.Value", contrast[[i]][1], "vs", contrast[[i]][2], sep = ".") # pval column
    logfcCol = paste("logFC", contrast[[i]][1], "vs", contrast[[i]][2], sep = ".") # logFC column
    
    logP = -log10(results[,pCol])
    logFC <- results[,logfcCol]
    sign = sign(logFC)
    
    if (rankMethod == "sign") { # geneList ranked by significance and FC direction
      
      metric = logP*sign
      geneList = metric
      names(geneList) = results$Geneid
      geneList = sort(geneList, decreasing = TRUE)
      
    } else if (rankMethod == "logFC") { # geneList ranked by significance and signal FC
      
      metric = logP*logFC
      geneList = metric
      names(geneList) = results$Geneid
      geneList = sort(geneList, decreasing=TRUE)
      
    }
    
    # GSEA
    set.seed(123)
    gsea[[i]] <- clusterProfiler::GSEA(geneList, TERM2GENE = gmt, verbose = FALSE, seed = T, minGSSize = minGSSize, maxGSSize = maxGSSize, pvalueCutoff = p.value, eps = 0, ...)
    
    
    # save excel with results
    headerStyle1 <- createStyle(halign = "center", valign = "center", textDecoration = "Bold",
                                wrapText = TRUE) 
    
    wb <- createWorkbook()
    for (j in 1:length(contrast[[i]])){ # Loop through positive and negative
      # Positive NES
      addWorksheet(wb, sheetName = paste("Enriched in", contrast[[i]][j]))
      if(j == 1){
        writeData(wb, j, gsea[[i]]@result[gsea[[i]]@result$NES > 0, -c(2,8)])
      } # POSITIVE NES
      if(j == 2){
        writeData(wb, j, gsea[[i]]@result[gsea[[i]]@result$NES < 0, -c(2,8)])
      } # NEGATIVE NES
      setColWidths(wb, sheet = j, cols = 1:9, widths = c(35,6,10,10,10,10,5,15,25))
      addStyle(wb, sheet = j, headerStyle1, rows = 1, cols = 1:9, gridExpand = TRUE)
      setRowHeights(wb, sheet = j, rows = 1, heights =30)
      
    }
    
    saveWorkbook(wb, file.path(gseaResDir, paste("GSEA", collectionName, contrast[[i]][1], "vs", contrast[[i]][2], "xlsx", sep=".")), overwrite = TRUE)
    
    message(paste0("\n --- Results for contrast ", contrast[[i]][1],".vs.", contrast[[i]][2]), " saved succesfully ---")
    
  }
  save(gsea, file = file.path(resultsDir, paste0("GSEA", collectionName, ".RData")))
  
  # make GSEA plots, calling the function
  if (plots) { 
    message("Drawing GSEA plots...")
    
    # check number of significant results for each contrast
    signif <- numeric(length = length(contrast))
    for (i in 1:length(contrast)){
      signif[i] <- length(which(gsea[[i]]@result$p.adjust < 0.05))
    }
    
    if (length(contrast) < 2 & length(which(signif > 0)) < 2) { # we need at least 2 contrasts with significant results to make the joinedDotplot
      JoinedDotplot = FALSE
      warning("Joined dotplot will be skipped: not enough contrasts")
      
    }
    
    # check which plots are to be made
    Barplot <- "Barplot" %in% plots2Make
    Dotplot <- "Dotplot" %in% plots2Make
    RunScoreSig <- "RunScoreSig" %in% plots2Make
    GCNetwork <- "GCNetwork" %in% plots2Make
    EnrichMAP <- "EnrichMAP" %in% plots2Make
    
    # call function to make plots
    makePlotsGSEA(gsea, contrast = contrast, collectionName = collectionName, resultsDir = resultsDir, plotTop = plotTop, plotTopRS = plotTopRS, p.adj = p.adj, JoinedDotplot = JoinedDotplot, Barplot = Barplot, Dotplot = Dotplot, RunScoreSig = RunScoreSig, GCNetwork = GCNetwork, EnrichMAP = EnrichMAP, contrast.order = contrast.order)
    
  }
  
  if (leadingEdges) {
    LEDir <- file.path(resultsDir,"leadingEdges")
    dir.create(LEDir, showWarnings = F)
    makeLeadingEdges(gsea, contrast = contrast, collectionName = collectionName, gseaDir = LEDir)
  }
  
  return(gsea)
  
  
}


##' Generates plots of Gene set enrichment analysis (GSEA)
##' 
##' Function that generates the plots to interpret GSEA results, some of them custom and some from the clusterProfiler package. All the following plots are generated by default: a general barplot, a dotplot, running scores, gene-concept networks, enrichment map, and a joined dotplot of common genesets in all conditions. 
##' 
##' @param gsea GSEA object generated with clusterProfiler. Returned from the `makeGSEA` function. 
##' @param contrast List with one vector of length 2 for each contrast.
##' @param collectionName Vector with the name of the collection.
##' @param resultsDir Character vector with output results directory. Default = working directory
##' @param Barplot Whether to perform the general barplot of significant results. If there are more than 30 significant gene sets up or down, only top 30 will be added. Default = TRUE
##' @param Dotplot Whether to perform GSEA dotplot. First 40 terms will be displayed. Default = TRUE
##' @param RunScoreSig Whether to perform the running scores of all significant terms. A maximum of 50 plots are generated by default. This number can be customized with parameter `plotTop`. Default = TRUE
##' @param GCNetwork Whether to perform gene-concept network. Displays top 5 gene sets by default. Default = TRUE
##' @param EnrichMAP Whether to perform enrichment map. It uses Jaccard index, and plots top 30 by default. Default = TRUE
##' @param GeneSetRS Name of gene set or list of gene sets to perform the running score instead of all significant. Default = NULL
##' @param JoinedDotplot Whether to perform the joined dotplot of common terms between all contrasts. Can also be obtained in separate function `makeJoinedDotplot`. Default = TRUE
##' @param plotTop Number of mamimum gene sets to be plotted in barplot and dotplot, and enrichment map. Default = 30
##' @param plotTopRS Number of maximum plots to perform the running score. Default = 50
##' @param p.adj Threshold of the adjusted p-value to be selected for plots. Default = 0.05
##' @param contrast.order Order in which the contrasts will be shown in the JoinedDotplot. If NULL, it will use the order as shown in contrast. Default = NULL
##' @param rankMethod Ranking method that was used for the GSEA. Either "sign", ranked by the log p-value and the sign of fold change, or "logFC", ranked by both log p-value and logFC. To be used for plot legend. Default = "sign"
##' 
##' @return Returns plots for each comparison in each correspondent subfolder at resultsDir 
##' @export
##' @import clusterProfiler
##' @import ggplot2
##' @import enrichplot
##' @import png
##' @import stringr
##' @import ggtangle


makePlotsGSEA <- function(gsea, contrast, collectionName = "", resultsDir = getwd(), Barplot = TRUE, Dotplot = TRUE, RunScoreSig = TRUE, GCNetwork = TRUE, EnrichMAP = TRUE, GeneSetRS = NULL, JoinedDotplot = TRUE, plotTop = 30, plotTopRS = 50, p.adj = 0.05, contrast.order = NULL, rankMethod = "sign", ...){
  
  contrast.t <- unlist(lapply(contrast, function(x) paste(x, collapse = ".vs.")))
  if(is.null(contrast.order)) {
    contrast.order = contrast.t
  } 
  
  for (i in 1:length(contrast)) {
    message("Plotting ", paste(contrast[[i]][1], "vs", contrast[[i]][2], sep="."))
    gseaResDir = file.path(resultsDir, paste("GSEA", collectionName, contrast[[i]][1], "vs", contrast[[i]][2], sep="."))
    if(!dir.exists(gseaResDir)) dir.create(gseaResDir)
    # dir.create(gseaResDir,showWarnings = F) #should be created already? make a ifnotexist or sth
    # Create a list with significant [[1]] positive and [[2]] negative results
    gsea.L = list()
    
    res = gsea[[i]]@result # results
    
    # positive NES
    gsea.L[[contrast[[i]][1]]] = gsea[[i]]
    gsea.L[[contrast[[i]][1]]]@result = res[res$NES > 0 & res$p.adjust < p.adj, ]
    # negative NES
    gsea.L[[contrast[[i]][2]]] = gsea[[i]]
    gsea.L[[contrast[[i]][2]]]@result = res[res$NES < 0 & res$p.adjust < p.adj, ]
    
    
    if (plotTop > 50) {
      warning("plotTop is set to a number higher than 50. This will difficult the interpretability and visualization of the plot. Please consider lowering the number.")
    }
    
    # Generate a general Barplot with NES and pval (maximum to plotTop up and down)
    up = gsea.L[[contrast[[i]][1]]]@result[, c("ID","NES","p.adjust")]
    up = up[1:ifelse(nrow(up) > plotTop, plotTop, nrow(up)),] # top 30 default (plotTOP)
    down = gsea.L[[contrast[[i]][2]]]@result[, c("ID", "NES", "p.adjust")]
    down = down[1:ifelse(nrow(down) > plotTop, plotTop, nrow(down)),]
    data = na.omit(rbind(up,down))
    # line wrapping for IDs that are too long
    data$ID = gsub("_", " ", data$ID)
    data$ID = str_wrap(data$ID, 70)
    # order by NES
    data$ID = factor(data$ID, data$ID[order(data$NES)], ordered = T)
    
    if (Barplot == TRUE) {
      if (nrow(data) != 0) {
        p = ggplot(data, aes(x = ID, y = NES, fill = p.adjust)) + 
          geom_bar(stat = "identity") +
          scale_fill_continuous(low = 'red', high = 'blue') +
          ylim(ifelse(min(data$NES) < (-2), min(data$NES), -2),
               ifelse(max(data$NES) > (2), max(data$NES), 2)) +
          coord_flip() +
          theme_bw() +
          ggtitle(paste0("GSEA of ", collectionName)) +
          theme(plot.title = element_text(hjust = 0.5, face="bold"), 
                axis.text.y = element_text(lineheight = 0.6, margin = margin(r = 5)))
        
        # check whether some lines are wrapped in order to increase size of the plot
        line_counts <- str_count(data$ID, "\n") + 1
        plot_height <- if (any(line_counts > 1)) {
          sum(line_counts) * 0.15
        } else {
          8
        }
        
        ggsave(file.path(gseaResDir, 
                         paste0("GSEA.", collectionName, ".BarplotNES.",
                                paste0(contrast[[i]][1], "vs", contrast[[i]][2]), ".png")),
               plot = p, width = 10, height = plot_height)
        
      } 
    }
    
    for (j in 1:length(contrast[[i]])) { # loop through positive and negative
      if (nrow(gsea.L[[j]]@result) != 0) { # If there are significant results; do plots
        
        # GSEA Dotplot
        if (Dotplot == TRUE) {
          
          p = clusterProfiler::dotplot(gsea.L[[j]], showCategory = plotTop, font.size = 8, label_format = 50, title =  paste0("Enriched in ", names(gsea.L)[j], "\n p.adjust < ", p.adj))
          
          ggsave(file.path(gseaResDir, paste0("GSEA.", collectionName, ".Top", plotTop, ".Dotplot.", names(gsea.L)[j], ".png")),
                 plot = p, width = 10, height = 10)
        }
        
        
        # GSEA Running score
        if (RunScoreSig == TRUE) {
          p = list()
          ## Create GSEA dir for Running Score results
          GSEA_score_dir = file.path(gseaResDir, "RunningScore")
          if(!dir.exists(GSEA_score_dir)) dir.create(GSEA_score_dir)
          
          ## Set max num. of plots to all DE or to plotTopRS
          n_plot = ifelse(nrow(gsea.L[[j]]@result) > plotTopRS, 
                          plotTopRS, 
                          nrow(gsea.L[[j]]@result))
          for (u in 1:n_plot){
            p = enrichplot::gseaplot2(gsea.L[[j]], color = "green", base_size = 4.5, 
                                      geneSetID = u, title = gsea.L[[j]]$Description[u])
            print(p)
            ggsave(filename = file.path(GSEA_score_dir, paste0(gsea.L[[j]]$Description[u], ".RunningScore.", names(gsea.L)[j],".png")),
                   width = 3,height = 2.5)
          }  
        }
        
        
        # GSEA Gene-Concept networks (top 5 gene sets by default)
        if (GCNetwork == TRUE) {
          # first without labels
          p = clusterProfiler::cnetplot(gsea.L[[j]], 
                                        foldChange = gsea.L[[j]]@geneList, 
                                        size_category = 1.2, 
                                        size_item = 1,
                                        size_edge = 0.3,
                                        layout = "kk",
                                        showCategory = 5, 
                                        node_label="none")
          
          # add gene labels
          p2 = p + geom_cnet_label(node_label = "item", size = 3, color = "gray40")
          # add category (geneset) labels in bold
          p3 = p2 + geom_cnet_label(node_label = "category", size=3.5, color="black", fontface = "bold")
          # add proper legend labeling, will depend on rank method
          LabelFC = ifelse(rankMethod == "logFC",
                           "-log(p.val)*log2FC",
                           "-log(p.val)*signFC") # by default should be else
          p3 = p3 + scale_color_gradient2(name = LabelFC, low = "blue", mid = "white", high = "red")
          
          # save plot
          ggsave(file.path(gseaResDir, paste0("GSEA.", collectionName, ".GeneConceptNetworks.",
                                              names(gsea.L)[j], ".png")), 
                 plot = p3, width = 10, height = 8, bg = "white")  
        }
        
        
        # GSEA EnrichmentMAP (Jaccard index. Plot top 30 by default (plotTop))
        if (EnrichMAP == TRUE) {
          if(nrow(gsea.L[[j]]@result) > 1) { # need results
            pt = enrichplot::pairwise_termsim(gsea.L[[j]], method = "JC", semData = NULL, showCategory = 200)
            
            if (nrow(pt@result) > 2) { # and at least 3 for plot
              p <- clusterProfiler::emapplot(pt, size_category = 0.5, showCategory = plotTop)
              # to reduce the size of labels like before
              label_layer <- which(sapply(p$layers, function(x) {
                "GeomTextRepel" %in% class(x$geom) || "GeomLabelRepel" %in% class(x$geom)
              }))
              # only modify if it exists (if empty fails)
              if (length(label_layer) > 0) {
                for (L in label_layer) {
                  p$layers[[L]]$aes_params$size <- 2.5
                }
              }
              ggsave(file.path(gseaResDir, paste0("GSEA.", collectionName, ".EnrichmentMAP.", 
                                                  names(gsea.L)[j], ".png")),
                     plot = p, width = 8, height = 6, bg = "white")
            } else {
              warning("Too few enriched terms for EnrichmentMap plot, skipping...")
            }
            
          }  
        }
        
      } else { # warn that there's no results 
          warning("No significant GSEA results for contrast: ", contrast.t[[i]], " (", names(gsea.L)[j], ")")
      
        }
      
      
      if (!is.null(GeneSetRS)) { # By default, significant running scores will be plotted. However, if a specific geneset is provided, will plot the running score plot of it (can be non-significant)
        if (length(GeneSetRS) == 1) { # if only one
          p = enrichplot::gseaplot2(gsea.L[[j]], color = "green", base_size = 4.5, 
                                    geneSetID = GeneSetRS, title = GeneSetRS)
          print(p)
          ggsave(filename = file.path(gseaResDir, paste0(GeneSetRS, ".RunningScore.", names(gsea.L)[j], ".png")),
                 width = 3, height = 2.5) 
        } else {
          for (gs in GeneSetRS) { # if more than one provided
            p = enrichplot::gseaplot2(gsea.L[[j]], color="green", base_size = 4.5, 
                                      geneSetID = gs, title = gs)
            print(p)
            ggsave(filename = file.path(gseaResDir, paste0(gs, ".RunningScore.", names(gsea.L)[j], ".png")),
                   width = 3, height = 2.5) 
          }
        }
        
      }
    }
  }
  
  if (JoinedDotplot == TRUE) {
    if (length(contrast) > 1) { # if contrast is 1 makes no sense to do the joined dotplot
      makeJoinedDotplot(gsea, contrast, contrast.order = contrast.order, p.adj = p.adj, collectionName = collectionName, resultsDir = resultsDir)
    } else {
      message("Joined dot plot aims to compare different contrasts. It will not be performed for just one contrast.")
    }
  }
  
  message("GSEA performed succesfully!")
  
  
}



#' Make GSEA dotplot
#' 
#' Generates a dotplot of GSEA results joining all comparisons when there are less than 100 significant genesets results, and an equivalent heatmap in Excel format for all cases 
#' 
#' @param gsea GSEA results object.
#' @param contrast List of vectors with each contrast to use
#' @param contrast.order Order in which the contrasts will be shown in plot. If NULL, it will use the order as shown in contrast. Default = NULL
#' @param p.adj Adjusted p-value threshold. Default = 0.05
#' @param collectionName Name of the collection used for the GSEA, e.g. "H" for Hallmark. Default = ""
#' @param resultsDir Output directory where results will be stored. Default = current working directory
#' @param sort_top If there are more than 100 significantly enriched genesets, only top 30 will be used for the plot. This parameter specifies whether the top genesets will be selected according to significance ("signif") or enrichment score ("NES"). Default = "signif"
#' @param wid Plot parameter: width of the device (in inches). If not specified, internal defaults will be used. Default = NULL
#' @param hei Plot parameter: height of the device (in inches). If not specified, internal defaults will be used. Default = NULL
#' @param excelFileName Name prefix to be used for the Excel heatmap file. Default = "GSEA.Results_HeatmapNES"

#' 
#' @import ggplot2
#' @import dplyr
#' @import openxlsx
#' @import stringr
#' @seealso [makeGSEA()]
#' 
#' @return Saves a dotplot of GSEA significant results in all conditions, as well as an excel with the same information
#' @export makeJoinedDotplot


makeJoinedDotplot <- function(gsea, contrast, contrast.order, p.adj = 0.05, collectionName = "", resultsDir = getwd(), sort_top = "signif", wid = NULL, hei = NULL, excelFileName = "GSEA.Results_HeatmapNES") {
  
  message("Creating joined dotplot...")
  contrast.t <- unlist(lapply(contrast, function(x) paste(x, collapse = ".vs.")))
  if(is.null(contrast.order)) {
    contrast.order = contrast.t
  } 
  
  result_list <- list()
  merged_df <- data.frame()
  # Loop to iterate each GSEA object
  for (i in 1:length(contrast)) {
    
    res_gsea <- gsea[[i]]@result
    result <- bind_rows(res_gsea[res_gsea$p.adjust < p.adj,], .id = "analysis")
    
    
    if (nrow(result) != 0) { # if there are significant results
      # add each result to results list
      result_list[[i]] <- result
      result_list[[i]]$analysis <- paste0(contrast[[i]][1],".vs.",contrast[[i]][2])
      df_heat <- result_list[[i]][,c("NES","p.adjust")]
      colnames(df_heat) <- paste0(colnames(df_heat), ".", result_list[[i]]$analysis[[1]])
      
      # create dataframe that will store all results for the later creation of "heatmap excel"
      if (nrow(merged_df) == 0) { # if empty it's the first added result
        merged_df <- df_heat
      } else { # if not the first, then we merge the existent with new one
        merged_df <- merge(merged_df, df_heat, by="row.names", all = TRUE)
        rownames(merged_df) <- merged_df$Row.names
        merged_df$Row.names <- NULL
      }
    }
    
  }
  
  # Add results into a single dataframe
  final_result <- bind_rows(result_list)
  final_result$analysis <- factor(final_result$analysis, levels = contrast.order) # show contrasts in order 
  
  # sort the columns of "heatmap excel" dataframe
  merged_df$Gene.Set.ID <- rownames(merged_df)
  res <- merged_df[, sort(colnames(merged_df))]
  
  
  # Remove the collection part from the GS name
  prefixes <- sapply(strsplit(final_result$ID, "_"), `[`, 1) # check all collection prefixes
  # remove the prefix only if all genesets are from hallmark, kegg or reactome
  if (length(unique(prefixes)) == 1) {
    if (unique(prefixes) %in% c("HALLMARK", "KEGG", "REACTOME")) {
      final_result$ID <- sub("HALLMARK_|KEGG_|REACTOME_", "", final_result$ID)
    }
  }
  
  # double line if the ID name is too long, needs to have space instead of _
  final_result$ID = gsub("_", " ", final_result$ID)
  final_result$ID = str_wrap(final_result$ID, 70)
  
  
  # Limit the number of significative results to make a plot.. and the letter size
  if (length(unique(final_result$ID)) <= 100) { # max 100 sig results to make the dotplot
    
    if (length(unique(final_result$ID)) > 50 & length(unique(final_result$ID)) <= 100) {
      charsize = 7
    } else {
      charsize = 10
    }
    
    if (is.null(wid) & is.null(hei)) { # set image size defaults
      hei = 8
      wid = ifelse(length(contrast) <= 4, 10, 12)
    }
    
    # make the Plot
    fplot <- final_result %>% 
      group_by(ID) %>% 
      ggplot(., aes(x = analysis, y = ID, size = p.adjust, col = NES))+ 
      geom_point()+ 
      scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', )+ 
      scale_size(range = c(1,5), trans = 'reverse')+ scale_y_discrete(limits = rev) +
      theme_bw() + 
      ggtitle(paste("Significant GSEA results in",collectionName))+
      theme(axis.title = element_blank(), 
            axis.text = element_text(size = charsize, face = 'bold'), 
            axis.text.x = element_text(angle = 45, hjust =1), 
            plot.title = element_text(size = 18, face = 'bold'), 
            strip.background = element_rect( color="black", fill="black", linewidth = 1.5, linetype="solid"), 
            strip.text = element_text(color = 'white', face = 'bold', size = 12),
            axis.text.y = element_text(lineheight = 0.6, margin = margin(r = 5)))
    
    
    ggsave(file.path(resultsDir, paste0("GSEA.", collectionName, ".Total.DotplotNES.png")),
           plot = fplot, width = wid, height = hei)
    
  } else { # more than 100 sig results  
    message("Too many significant results. Cannot perform a plot of all as it would be too confusing to interpret, therefore, only top 30 most ", ifelse(sort_top == "signif", "significant", "enriched"), " genesets of each comparison will be used. Please look at the equivalent excel file for comparing all enriched terms.")
    
    
    if (sort_top == "signif") {
      
      top30_df <- final_result %>%
        group_by(analysis) %>%
        slice_min(order_by = p.adjust, n = 30) %>%
        ungroup()
      
    } else if (sort_top == "NES") {
      
      top30_df <- final_result %>%
        group_by(analysis) %>%
        slice_max(order_by = abs(NES), n = 30) %>%
        ungroup()
      
    }
    
    if (is.null(wid) & is.null(hei)) { # set image size defaults
      hei = 10
      wid = ifelse(length(contrast) <= 4, 12, 14)
    }
    
    fplot30 <- top30_df %>% 
      group_by(ID) %>% 
      ggplot(., aes(x = analysis, y = ID, size = p.adjust, col = NES))+ 
      geom_point()+ 
      scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', )+ 
      scale_size(range = c(1,5), trans = 'reverse')+ scale_y_discrete(limits = rev) +
      theme_bw() + 
      ggtitle(paste("Significant GSEA results in", collectionName, "\ntop 30", ifelse(sort_top == "signif", "significant", "enriched"))) +
      theme(axis.title = element_blank(), 
            axis.text = element_text(size = 7, face = 'bold'), 
            axis.text.x = element_text(angle = 45, hjust =1), 
            plot.title = element_text(size = 18, face = 'bold'), 
            strip.background = element_rect( color="black", fill="black", linewidth = 1.5, linetype="solid"), 
            strip.text = element_text(color = 'white', face = 'bold', size = 12),
            axis.text.y = element_text(lineheight = 0.6, margin = margin(r = 5)))
    
    
    ggsave(file.path(resultsDir, paste0("GSEA.", collectionName, ".top30.", sort_top, ".DotplotNES.png")),
           plot = fplot30, width = wid, height = hei)
    
    
  }
  
  # Create Heatmap Excel with GSEA results for all conditions (NES + pvalue only)
  color.values <- res[,grep(colnames(res), pattern="NES", fixed = TRUE)]
  
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "GSEAresults")
  
  writeData(wb, 1, res[,-grep(colnames(res), pattern="p.adjust", fixed = TRUE)])
  
  headerStyle1 <- createStyle(fontSize = 10, textDecoration = "Bold",
                              wrapText = TRUE, textRotation = 90, halign = "center")
  
  addStyle(wb, sheet = 1, headerStyle1, rows = 1, cols = 1:length(colnames(res)),
           gridExpand = TRUE)
  
  number.col <- ncol(res) - dim(color.values)[2]
  
  # font size depending on pvalue levels
  pvalueLevel1 <- createStyle(fontSize = 8, numFmt = "0.00")
  pvalueLevel2 <- createStyle(fontSize = 10, numFmt = "0.00")
  pvalueLevel3 <- createStyle(fontSize = 12, numFmt = "0.00")
  pvalueLevel4 <- createStyle(fontSize = 14, numFmt = "0.00")
  
  for (i in 1:length(contrast)) {
    threshold_p <- paste0("p.adjust.", contrast.t[i])
    column_idx <- match(threshold_p, names(res))
    column_NES <- paste0("NES.", contrast.t[i])
    column_NES_idx <- match(column_NES, names(res))
    
    if (column_NES %in% names(res)) { # if there was no significant results for the comparison the column will not exist, so do it only if it does
      level1Rows <- which(res[,column_idx] < 0.05 & res[,column_idx] >= 0.005)
      level2Rows <- which(res[,column_idx] < 0.005 & res[,column_idx] >= 0.0005)
      level3Rows <- which(res[,column_idx] < 0.0005 & res[,column_idx] >= 0.00005)
      level4Rows <- which(res[,column_idx] < 0.00005)
      addStyle(wb, sheet = 1, style = pvalueLevel1, rows = level1Rows + 1, cols = column_NES_idx, gridExpand = TRUE)
      addStyle(wb, sheet = 1, pvalueLevel2, rows = level2Rows + 1, cols = column_NES_idx, gridExpand = TRUE)
      addStyle(wb, sheet = 1, pvalueLevel3, rows = level3Rows + 1, cols = column_NES_idx, gridExpand = TRUE)
      addStyle(wb, sheet = 1, pvalueLevel4, rows = level4Rows + 1, cols = column_NES_idx, gridExpand = TRUE)
      
    }
  }
  
  setRowHeights(wb, sheet = 1, rows = 1, heights = 150)
  setColWidths(wb, sheet = 1, cols = 1, widths = max(nchar(rownames(res))) + 1)
  setColWidths(wb, sheet = 1, cols = 2:number.col, widths = 8)
  
  freezePane(wb, sheet = 1 , firstRow = TRUE, firstCol = TRUE)
  
  # Heatmap :
  conditionalFormatting( wb,
                         sheet = 1,
                         cols = 2:as.numeric(dim(color.values)[2] + 1),
                         rows = 2:as.numeric(dim(color.values)[1] + 1),
                         rule = as.numeric(c(min(color.values, na.rm = TRUE), 0, max(color.values, na.rm = TRUE))),
                         style = c("blue","white", "red"),
                         type = "colorScale"
  )
  
  # legend of pvalue levels vs fontsize 
  addWorksheet(wb, sheetName = "Legend")
  
  legend.df <- data.frame(p.value = c(5e-2,5e-3,5e-4,5e-5),
                          fontSize = c("NES","NES","NES","NES"))
  
  options("openxlsx.borderColour" = "black")
  options("openxlsx.borderStyle" = "thin")
  writeData(wb, "Legend", legend.df, borderStyle = getOption("openxlsx.borderStyle", "thin"), startCol = 1, startRow = 1, borders = "columns")
  
  addStyle(wb, sheet = 2, pvalueLevel1, rows = 2, cols = 2, gridExpand = TRUE)
  addStyle(wb, sheet = 2, pvalueLevel2, rows = 3, cols = 2, gridExpand = TRUE)
  addStyle(wb, sheet = 2, pvalueLevel3, rows = 4, cols = 2, gridExpand = TRUE)
  addStyle(wb, sheet=2, pvalueLevel4, rows = 5, cols = 2, gridExpand = TRUE)
  
  
  saveWorkbook(wb, file.path(resultsDir, paste(excelFileName, collectionName, "xlsx", sep = ".")), overwrite = TRUE)
  
}

