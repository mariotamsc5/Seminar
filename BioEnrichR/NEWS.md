# BioEnrichR 1.0.0
This package is intended to execute an enrichment analysis on a RNA-seq dataset
## Functions or features
* `calculateDEGs()` calculates differentially expressed genes from gene expression count data. It utilizes the edgeR package for statistical analysis, filtering DEGs based on specified FDR and logFC thresholds. 
* `enrichGONE()` performs Gene Ontology enrichment analysis on a given list of genes. It reads DEGs data, maps gene symbols to ENTREZ IDs, and uses clusterProfiler to perform and visualize the GO enrichment analysis.
* `enrichKEG()` conducts KEGG pathway enrichment analysis. It processes input gene lists, performs KEGG pathway analysis, and generates various plots for interpretation.
## Documentation
The following files can be found with useful information about the package:
* DESCRIPTION: Package meta-data
* NAMESPACE: Functions and classes exported by the package
* LICENSE: License of the package
* NEWS (this document): Description of package revisions/versions

The following folders can be found with useful files and manuals of the package:
* R: R source code files
* Data: Test datasets
* man: Reference documentation
* Tests: Automated unit tests for package functions
* Vingettes: User manuals
