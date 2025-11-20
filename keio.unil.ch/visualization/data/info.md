# Information
keio.unil.ch provides an interactive platform for exploring FIA-TOF MS measurements of extracellular metabolite concentrations from batch cultures of E. coli grown in M9 glucose medium. The dataset covers 66 transporter-related knockout strains from the KEIO collection, alongside the wild-type reference strain (E. coli BW25113). Metabolites are annotated by m/z values using KEGG identifiers, names and metabolite classes.

The experiment was carried out in five independent batches, with six wild-type replicates in each batch. Each knockout strain was grown in three replicates within a single batch. Every data point represents the median standardized Z-score from four randomized FIA-TOF MS injections per sample. The cultivations were conducted by Snorre Sulheim in Sara Mitri's lab at UNIL and the mass-spec analyses were conducted by Peter F. Doubleday in Nicola Zamboni's lab at ETH. 

## Effect
We quantified the Effect of each gene knockout on a given metabolite as the ratio between the mean difference in extracellular concentration between the knockout and wild-type strains and the mean spread (mean absolute deviation) within each strain.

All notebooks, code, and processed data are available at
https://github.com/Mitri-lab/metabolite-release

Raw FIA-TOF MS data can be accessed on MassIVE under the dataset identifier MSV000097105.

For further details, see our preprint:
https://www.biorxiv.org/content/10.1101/2025.08.19.671024v1

If you find this resource useful, please share it—and cite our work when using these results.

## Contributors to keio.unil.ch
The online tool is hosted and maintained by Sara Mitri's lab at the Department of Fundamental Microbiology at the University of Lausanne (Switzerland).
The persons responsible for making this tool are:
- Snorre Sulheim (snorre.sulheim@unil.ch)
- Prajwal Padmanabha (prajwal.padmanabha@unil.ch)
- Eric Ulrich (eric.ulrich@unil.ch)






