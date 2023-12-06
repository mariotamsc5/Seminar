# This code is to draw the plots with enrichKEGG 

# Input_File_1: Path to the DEGs for the genes input in enrichKEGG
# Input_File_2: Path to the DEGs for the universe input in enrichKEGG
# Output_File: Path to the .xlsx file containing the results of enrichKEGG analysis.

enrichKEG<- function(Input_File_1,Input_File_2,Output_File) {
  
  
  library(clusterProfiler)
  library(edgeR)
  library(openxlsx)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggnewscale)
  library(ggupset)
  
  # Import file for gene attribute
  Df = read.xlsx(Input_File_1)
   
  # Importfile for universe attribute
  Df_Universe = read.xlsx(Input_File_2)
  
  Gene_List <- Df$logFC
  
  # Name the genes based on the symbol ID.
  names(Gene_List) <- Df[,1]
  
  # Sort the list in decreasing order
  Gene_List = sort(Gene_List, decreasing = TRUE)
  
  # Convert gene symbols to ENTREZID using bitr function
  Gene_List_Entrez<- bitr(names(Gene_List), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop=TRUE)
  
  # Check the converted gene list
  Gene_List_Entrez
  
  
  
  Universe_Gene_List = Df_Universe$logFC
  
  names(Universe_Gene_List) = Df_Universe[,1]
  
  Universe_Gene_List = sort(Universe_Gene_List,decreasing = TRUE)
  
  Universe_List_Entrez<- bitr(names(Universe_Gene_List), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop=T)
  
  Go_keggo <- enrichKEGG(gene = Gene_List_Entrez$ENTREZID,universe = universe_list_entrez$ENTREZID, organism  = "hsa", keyType = 'ncbi-geneid',pvalueCutoff = 0.05, qvalueCutoff = 0.1)
  
  
  
  # Bar plot of the data in the "Go_keggo" dataset
  show(barplot(Go_keggo))
  
  # Dt plot of the data in the "Go_keggo" dataset
  show(dotplot(Go_keggo))
  
  # Upset plot of the data in the "Go_keggo" dataset
  show(upsetplot(Go_keggo))
  
  # Heat map plot of the data in the "Go_keggo" dataset, with fold change values determined by the "gene_list" dataset
  show(heatplot(Go_keggo, foldChange = Gene_List))
  
  
  # Export the enrichKEGG results to an Excel file
  write.xlsx(Go_keggo, Output_File)
  
}

#Test
enrichKEG ("processed_data/DEGs_from_E-MTAB-2523.xlsx","processed_data/universe.xlsx","processed_data/kegging.xlsx")