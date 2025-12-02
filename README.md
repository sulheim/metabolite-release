![header](header.png)
# Summary
This repository contains data (except raw sequencing and metabolomics data), jupyter notebooks and other source code used in our recent project on microbial release rates. In this project we test the hypothesis that microbial release rates are negatively correlated with metabolites' inherent value. We also ask how important this factor is compared to other factors in explaining metabolite release rates and how microbial release rates change over a 4-6 week long evolution experiment with _E. coli_ kockout mutants.

# Repository folder structure
Detailed description of folder content is given in the README located in each folder. 

## Code
Matlab and python code used to parse and analyse data.

## Data
This folder contains all data (except sequencing data) used for the analyses and plots presented in the associated paper, including experimental data generated in this project and data collected from other sources. The raw sequencing data is available in the [NCBI SRA database](https://www.ncbi.nlm.nih.gov/sra) under the BioProject ID PRJNA1270783. The raw FIA-TOF mass-spectrometry data from the transporter KO experiment is available on [MassIVE](https://massive.ucsd.edu/ProteoSAFe/static/massive.jsp?redirect=auth) under the dataset identifier MSV000097105.

## keio.unil.ch
This folder contains the data and source code for the interactive web application hosted on [keio.unil.ch](https://keio.unil.ch). This web application can be used to explore the large dataset of untargeted exometabolome data from 66 _E. coli_ knockout mutants, each mutant devoid of a gene associated with metabolite transport across the inner or outer membrane. 

## Models
This folder contains the original and modified genome-scale metabolic models used in this project. This include the following species:
- _B. licheniformis_
- _C. glutamicum_
- _E. coli_
- _P. putida_
- _S. cerevisiae_

## Notebooks
All jupyter notebooks used for the analyses and plots presented in the associated paper.

# Reference
If you use these results, analyses or data, please cite this preprint:
[Sulheim, Snorre, Gunn Broli, Alisson Gillon, Julien S. Luneau, Andrew Quinn, Eric Ulrich, Margaret A. Vogel, Philipp Engel, and Sara Mitri. "Microbes release lower-value metabolites at higher rates." bioRxiv (2025): 2025-08.](https://www.biorxiv.org/content/10.1101/2025.08.19.671024v1)




