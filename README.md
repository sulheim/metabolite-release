![header](header.png)
# Summary
This repository contains all data, jupyter notebooks and other source code used in our recent project on microbial release rates. In this project we test the hypothesis that microbial release rates are shaped by the metabolites' inherent value. We also ask how important this factor is compared to other factors in explaining metabolite release rates and how microbial release rates change over a 4-6 week long evolution experiment with _E. coli_ kockout mutants.

# Repository folder structure
These are the key files and folders used in this project. 

## Data
- this_project
  - 1_e_coli_batch_cultures: Data from batch cultures of E. coli with galactose, L-malate or L-alanine as the main carbon source.
  - 2_keio_strains_screening: Data from cultivation and exometabolome screening of _E. coli_ KO mutants
  - 3_keio_strains_chemostat: Data from chemostat experiments with  _E. coli_ aceE and sucB KO mutants
  - 4_paired_metabolomics_live_dead: Paired intra- and extracellular metabolomics data and live/dead staining results
  - 5_div: Other data created in this project, e.g. files for mapping metabolite names to BiGG ids or to chemical properties
- paczia_2012: Original data from [Paczia et al. 2012](https://microbialcellfactories.biomedcentral.com/articles/10.1186/1475-2859-11-122) and files directly derived from this
- vila_2023: Original data from [Vila et al., 2023](https://doi.org/10.1101/2023.10.25.564019)

## Notebooks
- 1_release_rates
  - 1_sintef_2023: Estimate rates and metabolite values for _E. coli_ in galactose, L-malate and L-alanine
  - 2_paczia_et_al_2012: Estimate rates and metabolite values for _E. coli_, _B. licheniformis_, _C. glutamicum_, _S. cerevisiae_ in glucose
  - 3_vila_et_al: Estimate rates and metabolite values for _Enterobacter_ and _Pseudomonas_ environmental isolates
- 2_keio_strains_screening: Selection of KO strains and analyses of experimental data
- 9_other: Div notebooks used in this project

## Models
This folder contains the original and modified genome-scale metabolic models used in this project. This include the following species:
- _E. coli_
- _S. cerevisiae_
- _K. pneumoniae_
- _P. putida_
- _C. glutamicum_
- _B. subtilis_ / _B. licheniformis_
- _P. aeruginosa_


