#' Perform Gene Ontology Enrichment Analysis
#'
#' This function performs gene ontology (GO) enrichment analysis using the clusterProfiler package. It reads gene expression data from specified input files, converts gene symbols to ENTREZ IDs, performs enrichment analysis, generates various plots, and exports the results to an Excel file.
#'
#' @param Input_File_1 Path to the .xlsx file containing differentially expressed genes (DEGs) for GO enrichment analysis.
#' @param Input_File_2 Path to the .xlsx file containing universe genes for the analysis.
#' @param Output_File Path to save the results of the enrichGO analysis in .xlsx format.
#'
#' @return This function does not return a value but generates plots and writes the results of the enrichment analysis to an Excel file.
#' @import clusterProfiler edgeR openxlsx org.Hs.eg.db enrichplot ggnewscale ggupset
#' @export
#' @examples
#' enrichGONE("path/to/DEGs.xlsx", "path/to/universe_genes.xlsx", "path/to/enrichGONE_results.xlsx")
enrichGONE <- function(Input_File_1, Input_File_2, Output_File) {
  # ... [rest of the function code]
}
