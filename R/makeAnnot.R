#' Creates annotation file for a given genome
#' 
#' @param gtf Path to gtf file that was used for the alignment
#' @param genome Genome version that was used to create the count matrix. Example: mm10, GRCh38
#' @param counts Raw count matrix
#' @param annotDir Output directory where the annotations file will be stored
#' 
#' @return Saves an RData of the annotation file for the given genome version, at the provided output directory, with name "geneGTF.[genome].Annot.RData"
#' 
#' @import stringr
#' 
#' @export makeAnnot

makeAnnot <- function(gtf, genome, counts, annotDir) {
  
   if (genome %in% list.files(annotDir)) {
     stop(paste("The genome version",genome,"already has a previously created annotations file in the directory. You can read it using the function 'getAnnot'."))
     
   } else {
     message("Reading gtf file...")
     gtfData <- read.table(gtf, header = FALSE, sep = "\t")
     
     gtfSubset <- gtfData[gtfData$V3 == "exon",]  # gtf file has 'exon' and 'CDS' values
     
     # We want to keep:
     #   # gene id
     #   # initial position
     #   # ending position
     #   # strand (+, -)
     #   # chromosome
     #   # length
     
     geneList <- rownames(counts)
     N=length(geneList)
     iniList <- array(NA,N)
     endList <- array(NA,N)
     strandList <- array(NA,N)
     chrList <- array(NA,N)
     lenList <- array(NA,N)
     unfoundGenes <- list()
     
     message("Getting gene annotations from gtf file. This could take a while...")
     progress_bar = txtProgressBar(min=0, max=N, style = 3, char="=")
     
     for (i in 1:N) {
       gene <- geneList[i]
       geneLines <- gtfSubset[str_detect(gtfSubset$V9, paste("gene_id ", gene, ";", sep="")),]
       if (nrow(geneLines) > 0) {
         iniList[i] <- min(geneLines$V4)
         endList[i] <- max(geneLines$V5)
         strandList[i] <- geneLines$V7[1]
         chrList[i] <- geneLines$V1[1]
         lenList[i] <- endList[i]-iniList[i]
       } else {
         print(paste0(gene," gene not found"))
       }
       setTxtProgressBar(progress_bar, value = i)
     }
     close(progress_bar)
     geneListAnnot <- data.frame(Geneid = geneList,
                                 Chr = chrList,
                                 Start = iniList,
                                 End = endList,
                                 Strand = strandList,
                                 Length = lenList)
     rownames(geneListAnnot) <- geneList
     
     dir.create(file.path(annotDir,genome))
     #let's save it in the annotations for other projects
     output <- file.path(annotDir,genome,paste0("geneGTF.",genome,".Annot.RData"))
     save(geneListAnnot,file=output)
     message(paste("Annotations file for genome",genome,"was created successfully and saved at:",output))
   }
  
}