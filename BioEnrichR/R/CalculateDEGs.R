# This code is to calculate DEGs
# The DEGs are filtered by an FDR of <0.05 and an absolute logFC value of >1.3.

# Input_File: Path to the expression counts files.
# Output_File_1: Path to .xlsx file containing the list for most significant DEGs.
# Output_File_2: Path to .xlsx file containing the list of all DEGs unfiltered.
# Sample_Group: Vector to specifies the columns are control or disease.

cacuclateDEGs <- function(Input_File,Output_File_1,Output_File_2,Sample_Group) {
  library(edgeR)
  library(openxlsx)
  
  # Import counts.
  Gene_counts <- read.table(Input_File,header = T,as.is = T,row.names = 1,sep = "\t")
  
  dim(Gene_counts)
  head(Gene_counts)
  
  # Add more explanation.
  Sample_Groups <- factor(Sample_Group)
  
  # Create an edgeR list objet for DEGs.
  DEGs <- DEGList(counts = Gene_counts, group = Sample_Groups)
  
  # Filter low expression genes.
  # The genes with a minimum of 10 counts in at least 3 samples are retained as low-expression filtering.
  Filter <- filterByExpr(y = DEGs)
  DEG_Filter <- DEG[filter, keep.lib.sizes=F]
  
  # Caculating edgR's scaling factors for normalization
  DEG_Scal <- calNormFactors(object = DEG_Filter)
  
  # Fit the negative binomial model and estimate dispersion.
  DEG_Model <- estimateDisp(y = DEG_Scal)
  
  # Get P-values
  DEG_Stat <- exactTest(object = DEG_Model)
  
  # Caculate FDR
  DEG_Top = TopTags(object = DEG_Stat,n = "Inf")
  
  #Filter by FDR and logFC
  DEG_FDR <- DEG_Top[DEG_top$table$FDR<0.05,]
  
  DEG_FDR_logFC <- DEG_FDR[DEG_FDR$table$logFC<(-1.3)|DEG_FDR$table$logFC>1.3,]
  
  write.xlsx(DEG_FDR_logFC$table,Output_File_1,rowNames = TRUE)
  
  write.xlsx(DEG_Top$table,Output_File_2,rowNames = TRUE)
}

# Test
calculateDEGs("input_data/E-MTAB-2523.counts.txt","processed_data/DEGs_from_E-MTAB-2523.xlsx","processed_data/universe.xlsx",
              c(1,0,1,1,0,1,0,1,0,1,1,1,0,1,1,1,0,1))