#' Load annotation table
#' 
#' Gets annotation dataframe of a filtered counts matrix for an already created annotation file, that will be dependent on the genome species and version
#' 
#' @param counts Filtered count matrix.
#' @param genome Genome that was used to create the count matrix. Example: "mm10", "GRCh38"....
#' @param annotDir Directory where the annotation files are stored. The function will look for a subfolder with the genome name inside this directory.
#' 
#' @return Annotation dataframe of the genes included in the counts matrix.
#' 
#' @export getAnnot
#' @seealso [makeAnnot()]

getAnnot <- function(counts, genome, annotDir) {
  if (genome %in% list.files(annotDir)) {
    load(file.path(annotDir,genome,paste0("geneGTF.",genome,".Annot.RData")))
    dfAnnot <- geneListAnnot[rownames(counts),]
    return(dfAnnot)
  } else {
    stop(paste("The annotation for genome",genome,"could not be found in existent repertoire. Check if the path is correct, or create a new annotation file using the function 'makeAnnot'"))
  }
}
