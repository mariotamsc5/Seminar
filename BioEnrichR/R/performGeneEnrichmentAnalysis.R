#' Gene Enrichment Analysis
#'
#' @description
#' Prepares data from DEG analysis for Gene Enrichment Analysis using enrichGO from clusterprofiler.
#' It uses two files: one for DEGs (used in the 'genes' attribute of enrichGO) and one for the universe parameter, both including log fold change, p-adjustment, and FDR correction.
#' Generates barplot, dotplot, and upset plot for enrichment visualization.
#'
#' @param genesFilePath Path to the DEGs file for genes input in enrichGO.
#' @param universeFilePath Path to the DEGs file for universe input in enrichGO.
#' @param enrichmentResultsPath Path to the output .xlsx file containing enrichGO analysis results.
#'
#' @examples
#' library(DEGpathwayAnalysis)
#' performGeneEnrichmentAnalysis("processed_data/DEGs_from_E-MTAB-2523.xlsx",
#'                               "processed_data/universe.xlsx",
#'                               "processed_data/gene_enrichment_results.xlsx")
#'
#' @export
performGeneEnrichmentAnalysis <- function(genesFilePath, universeFilePath, enrichmentResultsPath) {

  library(clusterProfiler)
  library(edgeR)
  library(openxlsx)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggnewscale)
  library(ggupset)

  # Importing DEGs file for gene attribute
  genesDataFrame <- read.xlsx(genesFilePath)

  # Importing DEGs file for universe attribute
  universeDataFrame <- read.xlsx(universeFilePath)

  # Preparing gene list for enrichGO
  logFoldChangeList <- genesDataFrame$logFC
  names(logFoldChangeList) <- genesDataFrame[,1]
  logFoldChangeList <- sort(logFoldChangeList, decreasing = TRUE)
  entrezGeneList <- bitr(names(logFoldChangeList), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE)

  # Preparing universe gene list for enrichGO
  universeLogFoldChangeList <- universeDataFrame$logFC
  names(universeLogFoldChangeList) <- universeDataFrame[,1]
  universeLogFoldChangeList <- sort(universeLogFoldChangeList, decreasing = TRUE)
  entrezUniverseList <- bitr(names(universeLogFoldChangeList), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE)

  # Performing enrichGO analysis
  enrichGOAnalysis <- enrichGO(gene = entrezGeneList$ENTREZID, universe = entrezUniverseList$ENTREZID, OrgDb = "org.Hs.eg.db", keyType = 'ENTREZID', ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.1)

  # Generating and displaying plots
  show(barplot(enrichGOAnalysis))
  show(dotplot(enrichGOAnalysis))
  show(upsetplot(enrichGOAnalysis))
  show(heatplot(enrichGOAnalysis, foldChange = logFoldChangeList))

  # Exporting enrichGO results to an Excel file
  write.xlsx(enrichGOAnalysis, enrichmentResultsPath)
}

# Example usage
#performGeneEnrichmentAnalysis("processed_data/Significant_DEGs_E-MTAB-2523.xlsx", "processed_data/Complete_DEG_List.xlsx", "processed_data/gene_enrichment_results.xlsx")

