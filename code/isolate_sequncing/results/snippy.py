from os.path import join
import pandas as pd
import os

# Parsing annotated fixed mutations
dirs = os.listdir("../../data/isolate_sequencing/")
for dir in dirs:
    f = join("..", "..", "data", "isolate_sequencing", dir, "snippy", "snps.tab")
    df = pd.read_csv(f, sep="\t")
    df.to_csv(join("results", "fixed_mutations", dir + ".csv"), index=False)
