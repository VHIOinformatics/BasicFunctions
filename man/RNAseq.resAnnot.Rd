% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeRNAseqResults.R
\name{RNAseq.resAnnot}
\alias{RNAseq.resAnnot}
\title{Obtain results of RNAseq analysis}
\usage{
RNAseq.resAnnot(
  exprMat,
  annotMat,
  cond,
  fitMain = fit.main,
  contrast,
  species = "human",
  GO = TRUE,
  geneidCol = "Geneid",
  idType = "SYMBOL",
  resultsDir = getwd(),
  resAnnotFilename = "resultsAnnot",
  Excel = TRUE,
  excelFilename,
  pvalue = NULL,
  padj = 0.05,
  logFC = 1,
  add.colors = NULL
)
}
\arguments{
\item{exprMat}{Matrix of expression in CPM(E expression provided by voom), TMM or in counts, with rownames = Geneid.}

\item{annotMat}{Matrix with annotations obtained with `getAnnot` or `makeAnnot` functions.}

\item{cond}{Vector with sample conditions used in the contrasts, same order as colnames(exprMat).}

\item{fitMain}{Object obtained after contrasts.fit() function from limma package is applied. Default = fit.main}

\item{contrast}{List of vectors with each contrast to use.}

\item{species}{Species used to retrieve annotations (human or mouse). Default = "human"}

\item{GO}{Add annotation with GO. Default = TRUE}

\item{geneidCol}{Column of annotMat that contains gene identifier. Default = "Geneid"}

\item{idType}{type of gene identifier in column geneidCol. Possible values are "SYMBOL" or "ENSEMBL". Default = "SYMBOL"}

\item{resultsDir}{Output directory where results will be stored. Default = current working directory}

\item{resAnnotFilename}{Filename of the csv and rds files that will be saved, without the extension. Default = "resultsAnnot"}

\item{Excel}{Create an excel of the results. Default = TRUE}

\item{excelFilename}{Name of the excel file, without the extension}

\item{pvalue}{p-value threshold. Default = NULL, as it is assumed to use an adjusted p-value by default}

\item{padj}{Adjusted p-value threshold. Default = 0.05}

\item{logFC}{Abs(logFC) threshold. Default = 1}

\item{add.colors}{Colors to add if there are more than 20 contrasts. Default = NULL}
}
\value{
Results of RNA-seq analysis. Returns results object, and saves outputs in csv and rds format to resultsDir. By default, results are also saved in Excel format.
}
\description{
Creates results of differential expression RNAseq analysis, which may including GO annotations and an excel with the total and per comparison results
}
\seealso{
[makeExcelResults(), getAnnot(), makeAnnot()]
}
