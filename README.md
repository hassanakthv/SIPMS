![SIP-MS Logo](https://github.com/hassanakthv/SIPMS/assets/43888767/70437bd0-88f8-4591-8b08-c4f5215e6713)
# Species Identification and Prediction by Mass Spectrometry (SIP-MS) 

## Species Identification and Prediction Method

SIP-MS leverages shotgun proteomics techniques to offer collagenous peptide-based species identification. It stands on two foundational pillars: a machine learning method classifier (a Random Forest classifier) with species-specific peptide sequences and abundances, and a correlation classifier that considers all informative peptides in a dataset.
***

## Species Search Engine (SSE)

SIP-MS is integrated as the back-end algorithm for the current SSE GUI, which, as of January 2024, encompasses 8 species. SSE takes into account various factors provided by SIP-MS results. It then either yields a prediction score for the submitted sample to identify the species or offers a similarity score for a sample in cases where SSE determines that the species for the submitted samples is not in the current database.


Here you can see the current species inventory
![SSE_Description of Current Database](https://github.com/hassanakthv/SIPMS/assets/43888767/b38933a0-56c3-4b79-b6b5-5944f864477b)
***
## Performing Prediction on a local server
Currently, SSE can be utilized locally with the following command in R versions higher than 4.0:

```
install.packages("remotes")
install.packages("shiny")

library(shiny)
library(remotes)

install_github("hassanakthv/SIPMS")
runGitHub(repo = "hassanakthv/SIPMS",subdir = "R")
```
The peptide file should have a column named "Peptide" as the first column and the rest of the columns should denote the samples names and contain the peptide abundances.
