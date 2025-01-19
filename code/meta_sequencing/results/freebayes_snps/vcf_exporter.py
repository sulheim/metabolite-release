import vcfpy
from os.path import join, exists
from os import listdir
import pandas as pd

"""This scrpts parses the vcf file created from freebayes. It showed that the anlaysis is identical to snippy.
This is a good test that shows that the mapping, freebayes calling and parsing is solid which is 
what I will use for the metagenomic data."""
dfs = []
for d in listdir(
    "/work/FAC/FBM/DMF/smitri/evomicrocomm/seq_snorre/data/meta_sequencing"
):
    f = join(
        "/work/FAC/FBM/DMF/smitri/evomicrocomm/seq_snorre/data/meta_sequencing",
        d,
        "var.annotated.vcf",
    )
    if exists(f):
        reader = vcfpy.Reader.from_path(f)
        for record in reader:
            # Iterating over each variant in case different variants
            # are present at the same position
            for i, r in enumerate(record.ALT):
                snp = {}
                # Chromosome
                snp["chrom"] = record.CHROM
                # Position
                snp["pos"] = record.POS
                # Variant quality
                snp["qual"] = record.QUAL
                # Sequencing depth
                snp["depth"] = record.INFO["DP"]
                # Variant frequency
                snp["freq"] = record.INFO["AO"][i] / record.INFO["DP"]
                # Alt sequence
                snp["alt"] = r
                # Alt depth
                snp["alt_count"] = record.INFO["AO"][i]
                # Reference sequence
                snp["ref"] = record.REF
                # Mutation type
                snp["type"] = record.INFO["TYPE"][i]
                # Mutation length for INS and DEL
                snp["len"] = record.INFO["LEN"][i]
                # Coding or non-coding
                if "ANN" in record.INFO.keys():
                    snp["eff"] = record.INFO["ANN"][0].split("|")[1]
                    snp["gene"] = record.INFO["ANN"][0].split("|")[4]
                key = ".".join(
                    [
                        snp["chrom"],
                        str(snp["pos"]),
                        str(d),
                    ]
                )
                snp["linegroup"] = key
                # Quality cutoff
                df = pd.concat(
                    [
                        pd.DataFrame(snp, index=[0]),
                        pd.DataFrame(data={"sample": [d]}, index=[0]),
                    ],
                    axis=1,
                )
                dfs.append(df)

out = pd.concat(dfs)
for sample in set(out["sample"]):
    df = out[out["sample"] == sample]
    df.to_csv(sample + ".csv", index=False)
    filter = (df["freq"] >= 0.05) & (df["qual"] >= 1000) & (df["alt_count"] >= 5)
    df = df[filter]
    df.to_csv(sample + ".filtered.csv", index=False)

filter = (out["freq"] >= 0.05) & (out["qual"] >= 1000) & (out["alt_count"] >= 5)
out = out[filter]
out.to_csv("all_samples.filtered.csv", index=False)
