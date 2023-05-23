##' plot ORA results, internal function
##' 
##' @param enrichment results of an ORA enrichment performed with function makeORA
##' @param plotName the name to be given to that plot
##' @param collectionName name of the mSigDB collection. 
##' @param resultsDir Character vector with output results directory
##' @param plotTop Number of maximal gene sets to be plotted in barplot and dotplot. Defalult = 50
##' @param plotP.adjust threshold of the adjusted p-value to be selected for plots. Default = 0.05

##' @import clusterProfiler
##' @import enrichplot
##' @import gridExtra
##' @import png
##' @import ggplot2

plotORA <- function(enrichment, plotName, collectionName = "", resultsDir = getwd(), plotTop = 50, plotP.adjust = 0.05) {
  
  
  enrichment@result = enrichment@result[enrichment@result$p.adjust < p.adjust,]
  fullName <- paste(collectionName, plotName, sep = ".")
  
  if (nrow(enrichment@result) > 0) { # If there are significant results; do plots
    
    enrichment@result$Description=strtrim(enrichment@result$Description, 80) # maximum label length
    p = barplot(enrichment, showCategory=plotTop, font.size = 6, title = paste0("Top ", plotTop, " ", collectionName, " enriched in ", plotName, " , p.adjust < ", p.adjust))
    ggsave(file.path(resultsDir, paste("Enrichment", fullName, "Barplot", "png", sep = ".")), plot = p, width = 9, height = ifelse(nrow(enrichment@result) > 5, 8, 4))
    
    p = clusterProfiler::dotplot(enrichment, showCategory = plotTop, font.size = 6, title = paste0("Top ", plotTop, " ", collectionName, " enriched in ", plotName, " , p.adjust < ", p.adjust))
    ggsave(file.path(resultsDir, paste("Enrichment", fullName, "Dotplot", "png", sep = ".")), plot = p, width = 9, height = ifelse(nrow(enrichment@result) > 5, 8, 4))
    
    p = clusterProfiler::cnetplot(enrichment, cex.params = list(category_node = 0.7, category_label = 0.7, gene_label = 0.5), layout = "kk", showCategory = 10)
    ggsave(file.path(resultsDir, paste("Enrichment", fullName, "GeneConceptNetworks", "png", sep = ".")), plot = p, width = 9, height = 8)
    
    if (nrow(enrichment@result) > 1) {
      
      pt = enrichplot::pairwise_termsim(enrichment, method = "JC", semData = NULL, showCategory = 200)
      
      p <- clusterProfiler::emapplot(pt, cex.params = list(category_label = 0.5), showCategory = 30)
      ggsave(file.path(resultsDir, paste("Enrichment", fullName, "EnrichmentMAP", "png", sep = ".")), plot = p, width = 9, height = ifelse(nrow(enrichment@result) > 5, 8, 4))
    }
  }
  
}