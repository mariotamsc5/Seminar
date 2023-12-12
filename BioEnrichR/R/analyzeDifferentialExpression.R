#' Analyze Differential Gene Expression
#'
#' @description
#'
#' This function filters and normalizes RNA-seq expression data, calculates
#' statistically significant differentially expressed genes (DEGs), and exports
#' them to an Excel file. It filters DEGs based on FDR < 0.05 and absolute logFC value > 1.3.
#' The low expression filtering parameters are set to keep genes with at least 10 counts
#' in at least 3 samples.
#'
#' @param expressionDataPath Path to the expression counts file.
#' @param significantDEGsOutputPath Path to the .xlsx file that will be created,
#'        containing a list of most significant DEGs.
#' @param allDEGsOutputPath Path to the .xlsx file used as a universe in the
#'        other functions in the program, containing a list of all DEGs unfiltered.
#' @param groupClassification A vector of ones and zeros specifying which samples
#'        (columns) are control/healthy (0) and which are disease (1).
#' @import edgeR
#' @import openxlsx
#' @examples
#' analyzeDifferentialExpression("E-MTAB-2523.counts.txt",
#'                               "Significant_DEGs_E-MTAB-2523.xlsx",
#'                               "Complete_DEG_List.xlsx",
#'                               c(1,0,1,1,0,1,0,1,0,1,1,1,0,1,1,1,0,1))
#' @export


analyzeDifferentialExpression <- function(expressionDataPath, significantDEGsOutputPath, allDEGsOutputPath, groupClassification) {

  library(edgeR)
  library(openxlsx)

  # Reading gene expression counts
  expressionCounts <- read.table(expressionDataPath, header=TRUE, as.is=TRUE, row.names=1, sep="\t")

  # Group classification provided by the user
  groupFactors <- factor(groupClassification)

  # Creating an edgeR DGEList object
  differentialExpressionList <- DGEList(counts = expressionCounts, group = groupFactors)

  # Filtering lowly expressed genes
  expressionFilter <- filterByExpr(y = differentialExpressionList)
  filteredExpressionList <- differentialExpressionList[expressionFilter, keep.lib.sizes=FALSE]

  # Normalizing using edgeR's scaling factors
  normalizedExpressionList <- calcNormFactors(object = filteredExpressionList)

  # Modeling and estimating dispersion
  modeledExpressionList <- estimateDisp(y = normalizedExpressionList)

  # Conducting exact test for p-values
  expressionStatistics <- exactTest(object = modeledExpressionList)

  # Identifying top tags (DEGs) and applying FDR and logFC thresholds
  topExpressionTags <- topTags(object = expressionStatistics, n = "Inf")
  filteredDEGsByFDR <- topExpressionTags[topExpressionTags$table$FDR < 0.05,]
  significantDEGs <- filteredDEGsByFDR[filteredDEGsByFDR$table$logFC < (-1.3) | filteredDEGsByFDR$table$logFC > 1.3,]

  # Exporting results to Excel files
  write.xlsx(significantDEGs$table, significantDEGsOutputPath, rowNames = TRUE)
  write.xlsx(topExpressionTags$table, allDEGsOutputPath, rowNames = TRUE)
}

