#' KEGG Pathway Enrichment Analysis
#'
#' @description
#' Conducts KEGG Pathway Enrichment Analysis on DEGs using enrichKEGG from the clusterProfiler package.
#' Utilizes filtered and log fold change adjusted DEG data files for gene input and universe input in enrichKEGG.
#' Generates barplot, dotplot, and upset plot for enrichment visualization.
#'
#' @param genesDataPath Path to the DEGs file for gene input in enrichKEGG.
#' @param universeDataPath Path to the DEGs file for universe input in enrichKEGG.
#' @param keggResultsPath Path to the .xlsx file that will be created, containing enrichKEGG analysis results.
#'
#' @examples
#' library(DEGpathwayAnalysis)
#' performKEGGAnalysis("processed_data/DEGs_from_E-MTAB-2523.xlsx",
#'                     "processed_data/universe.xlsx",
#'                     "processed_data/kegg_analysis_results.xlsx")
#'
#' @export
performKEGGAnalysis <- function(genesDataPath, universeDataPath, keggResultsPath) {

  library(clusterProfiler)
  library(edgeR)
  library(openxlsx)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggnewscale)
  library(ggupset)

  # Importing DEGs data for genes
  genesDataFrame <- read.xlsx(genesDataPath)

  # Importing DEGs data for universe
  universeDataFrame <- read.xlsx(universeDataPath)

  # Preparing gene list for enrichKEGG
  genesLogFC <- genesDataFrame$logFC
  names(genesLogFC) <- genesDataFrame[,1]
  genesLogFC <- sort(genesLogFC, decreasing = TRUE)
  genesEntrezIDs <- bitr(names(genesLogFC), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE)

  # Preparing universe gene list for enrichKEGG
  universeLogFC <- universeDataFrame$logFC
  names(universeLogFC) <- universeDataFrame[,1]
  universeLogFC <- sort(universeLogFC, decreasing = TRUE)
  universeEntrezIDs <- bitr(names(universeLogFC), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE)

  # Conducting enrichKEGG analysis
  keggAnalysisResults <- enrichKEGG(gene = genesEntrezIDs$ENTREZID, universe = universeEntrezIDs$ENTREZID, organism  = "hsa", keyType = 'ncbi-geneid', pvalueCutoff = 0.05, qvalueCutoff = 0.1)

  # Plotting results
  show(barplot(keggAnalysisResults))
  show(dotplot(keggAnalysisResults))
  show(upsetplot(keggAnalysisResults))
  show(heatplot(keggAnalysisResults, foldChange = genesLogFC))

  # Exporting results to Excel file
  write.xlsx(keggAnalysisResults, keggResultsPath)
}

# Example usage
#performKEGGAnalysis("processed_data/Significant_DEGs_E-MTAB-2523.xlsx", "processed_data/Complete_DEG_List.xlsx", "processed_data//kegg_analysis_results.xlsx")
