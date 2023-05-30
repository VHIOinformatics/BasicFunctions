# BasicFunctions
Basic functions we use on our projects for bioinformatics data analysis.

## Installation

```r
library(devtools)
install_github("VHIOinformatics/BasicFunctions")
```

## Package contents

* **createTargets.R**: Creates a file of a table of targets of the project. Includes information such as the sample names, treatment, status, and condition (the column information will vary depending on the project). The function `createTargets` saves a txt file of the generated table.

* **normTMM.R**: The function `normTMM` returns a count matrix with TMM normalized values.

* **rawMatrix2Excel.R**: The function `rawMatrix2Excel` converts a count matrix to Excel format, ready to send to investigators who may ask for it.   

* **getAnnot.R**: The function `getAnnot` loads the dataframe of an existing annotation file, specific of the genome version. It will check if a file with the structure "geneGTF.<genome>.Annot.RData" is present at provided directory, and return an annotation table of genes present in the count matrix provided as input (can be filtered or unfiltered).

* **makeAnnot.R**: If the annotation file has not been created for an specific genome, the function `makeAnnot` will generate it, and save it in a folder with the genome version name, and with the structure "geneGTF.<genome>.Annot.RData". The input count matrix should be raw (not filtered) to ensure this annotation file can be used in other projects.

* **makeRNAseqResults.R**: Script to obtain main results of RNAseq analysis. Includes functions `RNAseq.resAnnot`, which creates a table of differential expression results and saves it in "rds" and "csv" formats, with or without GO annotations, and `makeExcelResults`, which creates an Excel file with the results when called by the `RNaseq.resAnnot` or independently.

* **makeCls.R**: Creates the vector of conditions for GSEA in Cls format.

* **makeGct.R**: Creates gct file from a matrix with the format needed for GSEA.

* **makeRnk.R**: Creates rnk file from p-values and logFC to use in GSEA.

* **makeGSEA.R**: Includes functions `makeGSEA`, `makePlotsGSEA` and `makeJoinedDotplot`. `makeGSEA` performs the Gene Set Enrichment Analysis from the results object generated with `RNAseq.resAnnot`. It saves an Excel with significative genes in each gene set for each comparison, and calls the `makePlotsGSEA` function to draw the plots (Barplot, Dotplot, RunningScores, Gene-concept networks, Enrichment map). `makePlotsGSEA` additionally calls for `makeJoinedDotplot`, which may also be called independently, to generate a dotplot of all conditions combined (in Excel format and as an image if there are less than 100 significant gene set results).

* **makeORA.R**: Includes functions `makeORA`, which runs the Over Representation Analysis, and `plotORA`, which is called by the `makeORA` function to generate the plots of the results.
