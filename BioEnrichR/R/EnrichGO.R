# This code is to draw the plots with enrichGO 

# Input_File_1: Path to the DEGs for the genes input in enrichGO.
# Input_File_2: Path to the DEGs for the universe input in enrichGO.
# Output_File: Path to the .xlsx file containing the results of enrichGO analysis.

enrichGO <- function(Input_File_1,Input_File_2,Output_File) {
  
  
  library(clusterProfiler)
  library(edgeR)
  library(openxlsx)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggnewscale)
  library(ggupset)
  
  # Import file for the gene attribute.
  Df = read.xlsx(Input_File_1)
  
  # Import file for universe attribute
  Df_Universe = read.xlsx(Input_File_2)
  
  Gene_List <- Df$logFC
  
  
  # Name the genes based on the symbol ID.
  names(Gene_List) <- Df[,1]
  
  # Sort the list in decreasing order.
  Gene_List = sort(Gene_List, decreasing = TRUE)
  
  # Convert gene symbols to ENTREZID using bitr function
  Gene_List_Entrez<- bitr(names(Gene_List), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db",drop=TRUE)
  
  # check the converted gene list
  Gene_List_Entrez

  
  
  Universe_Gene_List = Df_Universe$logFC

  names(Universe_Gene_List) = Df_Universe[,1]
  
  Universe_Gene_List = sort(Universe_Gene_List,decreasing = TRUE)
  
  Universe_List_Entrez <- bitr(names(Universe_Gene_List), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop=T)
  
  
  Go_Enrich <- enrichGO(gene = Gene_List_Entrez$ENTREZID,universe = Universe_List_Entrez$ENTREZID, OrgDb = "org.Hs.eg.db" , keyType = 'ENTREZID',ont = "BP",pvalueCutoff = 0.05, qvalueCutoff = 0.1)
  
  
  
  # Bar plot for the data in the "Go_Enrich" dataset
  show(barplot(Go_enrich))
  
  # Dot plot for the data in the "Go_Enrich" dataset
  show(dotplot(Go_enrich))
  
  # Upset plot of the data in the "Go_Enrich" dataset
  show(upsetplot(Go_enrich))
  
  # Heat map plot of the data in the "Go_Enrich" dataset
  show(heatplot(Go_enrich, foldChange = Gene_List))
  
  # Export the enrichGO results to the Excel file.
  write.xlsx(Go_enrich, output_file)
  
}

#Test
enrichGO ("processed_data/DEGs_from_E-MTAB-2523.xlsx","processed_data/universe.xlsx","processed_data/enrichGONE_results.xlsx")