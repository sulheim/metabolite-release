# keio.unil.ch
[keio.unil.ch](https://keio.unil.ch) provides an interactive platform for exploring FIA-TOF MS measurements of extracellular metabolite concentrations from batch cultures of E. coli grown in M9 glucose medium. The dataset covers 66 transporter-related knockout strains from the KEIO collection, alongside the wild-type reference strain (E. coli BW25113). Metabolites are annotated by m/z values using KEGG identifiers, names and metabolite classes.

The experiment was carried out in five independent batches, with six wild-type replicates in each batch. Each knockout strain was grown in three replicates within a single batch. Every data point represents the median standardized Z-score from four randomized FIA-TOF MS injections per sample. The cultivations were conducted by Snorre Sulheim in Sara Mitri's lab at UNIL and the mass-spec analyses were conducted by Peter F. Doubleday in Nicola Zamboni's lab at ETH.

## Run it at home
The application is developed using Dash in python, and can easily be run locally: download this folder, navigate to `keio.unil.ch/visualization` and start the application by invoking the command:`python dashboard.py`
