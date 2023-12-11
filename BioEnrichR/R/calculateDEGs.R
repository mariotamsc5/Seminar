#' Calculate Differentially Expressed Genes (DEGs)
#'
#' This function calculates differentially expressed genes (DEGs) based on expression count data.
#' DEGs are identified using the edgeR package and filtered based on False Discovery Rate (FDR)
#' and log fold change (logFC) thresholds. The results are written to Excel files.
#'
#' @param Input_File A string specifying the path to the expression counts file.
#' @param Output_File_1 A string specifying the path to the .xlsx file where the most significant DEGs will be saved.
#' @param Output_File_2 A string specifying the path to the .xlsx file where all DEGs will be saved, unfiltered.
#' @param Sample_Group A numeric vector specifying which columns in the count data are control or disease.
#' @param FDR_Threshold A numeric value for the FDR cutoff for filtering DEGs. Default is 0.05.
#' @param logFC_Threshold A numeric value for the log fold change cutoff for filtering DEGs. Default is 1.3.
#' @param verbose A boolean indicating whether to print additional information (dimensions and head of gene counts). Default is TRUE.
#'
#' @return A list containing two elements: Significant_DEGs (DEGs filtered based on FDR and logFC) and All_DEGs (all DEGs unfiltered).
#' @export
#' @examples
#' result <- calculateDEGs("path/to/expression_counts.txt",
#'                         "path/to/significant_DEGs.xlsx",
#'                         "path/to/all_DEGs.xlsx",
#'                         c(1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1),
#'                         verbose = FALSE)
calculateDEGs <- function(Input_File, Output_File_1, Output_File_2, Sample_Group, FDR_Threshold = 0.05, logFC_Threshold = 1.3, verbose = TRUE) {
  # ... [rest of the function code]
}

