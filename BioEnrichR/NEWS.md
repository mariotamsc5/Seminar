# BioEnrichR 1.0.0
This package is intended to execute an enrichment analysis on a RNA-seq dataset
## Functions or features
* `calculateDEGs()` calculate DEGs filtered by an FDR of <0.05 and an absolute logFC value of >1.3. The output given is the path to .xlsx file containing the list for most significant DEGs and of all DEGs unfiltered.
* `enrichGONE()` executes the plots with enrichGO. The output given is the path to the .xlsx file containing the results of enrichGO analysis.
* `enrichKEG()` executes the plots with enrichKEG. The output given is the path to the .xlsx file containing the results of enrichKEGG analysis.
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
