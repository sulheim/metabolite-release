# Information
keio.unil.ch is an online tool for browsing FIA-TOF MS data on extracellular metabolite concentrations across batch cultivations of different _E. coli_ strains in M9 glucose medium. The dataset comprises 66 knockout (KO) strains from the KEIO collection and the KEIO collection wild-type (WT) reference strain (_E. coli_ BW25113). Each KO strain lacks one gene associated with metabolite transport. Metabolites are annotated based on m/z values to KEGG annotations, metabolite names and metabolite classes. Due to the size of the experiment, it was conducted in five seperate batches and the WT was cultivated is six replicates in each batch. All KO strains was cultivated in three replicates in one of the five batches. Each datapoint shown here is the median standardized Z-score of 4 FIA-TOF MS injections of each sample. Sample injections were randomized within each plate.

## Effect
We quantify the _Effect_ of a gene knockout on a specific metabolite as the ratio of the mean difference in extracellular concentration between the corresponding KO strain and the WT strain to the average absolute variability in extracellular concentrations within each strain. 

Notebooks, code snippets and data used in this work are found at https://github.com/Mitri-lab/metabolite-release. The raw FIA-TOF MS data is available on MassIVE under then dataset identifier MSV000097105.

For more details we refer to the preprint: https://www.biorxiv.org/content/10.1101/2025.08.19.671024v1. If you like this work please share it. If you use these results please cite us.

## Contributors to keio.unil.ch
The online tool is hosted and maintained by Sara Mitri's lab at the Department of Fundamental Microbiology at the University of Lausanne (Switzerland).
The persons responsible for making this tool are:
- Snorre Sulheim (snorre.sulheim@unil.ch)
- Prajwal Padmanabha (prajwal.padmanabha@unil.ch)
- Eric Ulrich (eric.ulrich@unil.ch)






