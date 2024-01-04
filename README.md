![SIP-MS Logo](https://github.com/hassanakthv/SIPMS/assets/43888767/70437bd0-88f8-4591-8b08-c4f5215e6713)
# Species Identification and Prediction by Mass Spectrometry (SIP-MS) 

*Species Identification and Prediction Method*

SIP-MS uses shotgun proteomics techniques to provide collagenous peptide-based species identification. SIP-MS is built upon two pillars: a machine learning method classifier (a Random Forest classifier) with species-specific peptide sequences and abundances and a correlation classifier that considers all the informative peptides in a dataset. 
***

*Species Search Engine (SSE)*

SIP-MS is embedded as the back-end algorithm for the current SSE GUI which as of Jan 2024 comprises 8 species. SSE takes into account several factors that are given as the results in SIP-MS and then either gives a prediction score for the submitted sample to identify the species or provides a similarity score for a sample (in case that the verdict from SSE is that the species for the submitted samples is not in the current database).


Here you can see the current species inventory
![SSE_Description of Current Database](https://github.com/hassanakthv/SIPMS/assets/43888767/b38933a0-56c3-4b79-b6b5-5944f864477b)
***
# Performing Prediction on a local server
Currently, SSE can be used locally using the below command in the R version higher than 4.0

``` 
__shiny::runGitHub(repo = "hassanakthv/SIPMS",subdir = "R")__
``` 
