#' Perform KEGG Pathway Enrichment Analysis
#'
#' This function performs Kyoto Encyclopedia of Genes and Genomes (KEGG) pathway enrichment analysis on a set of genes. It utilizes the clusterProfiler package to analyze gene lists against KEGG pathways, generating various plots and exporting the results to an Excel file. The function reads gene expression data, converts gene symbols to ENTREZ IDs, and performs the enrichment analysis.
#'
#' @param Input_File_1 Path to the .xlsx file containing differentially expressed genes (DEGs) for KEGG enrichment analysis.
#' @param Input_File_2 Path to the .xlsx file containing universe genes for the analysis.
#' @param Output_File Path to save the results of the enrichKEGG analysis in .xlsx format.
#'
#' @return This function does not return a value but generates plots and writes the results of the enrichment analysis to an Excel file.
#' @import clusterProfiler edgeR openxlsx org.Hs.eg.db enrichplot ggnewscale ggupset
#' @export
#' @examples
#' enrichKEG("path/to/DEGs.xlsx", "path/to/universe_genes.xlsx", "path/to/enrichKEGG_results.xlsx")
enrichKEG <- function(Input_File_1, Input_File_2, Output_File) {
  # ... [rest of the function code]
}
