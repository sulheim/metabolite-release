import vcfpy
from os.path import join, exists
from os import listdir
import pandas as pd
from BCBio import GFF
import numpy as np

"""This scrpts parses the vcf file created from freebayes. It showed that the anlaysis is identical to snippy.
This is a good test that shows that the mapping, freebayes calling and parsing is solid which is 
what I will use for the metagenomic data."""


def parse_gff():
    gff_f = "/work/FAC/FBM/DMF/smitri/evomicrocomm/seq_snorre/data/references_sequencing/reference.gff"
    gene_dict = {}
    with open(gff_f, "r") as handle:
        for record in GFF.parse(handle):
            for feature in record.features:
                if "gene" in feature.qualifiers.keys():
                    gene = feature.qualifiers["gene"][0]
                else:
                    gene = ""
                if "product" in feature.qualifiers.keys():
                    product = feature.qualifiers["product"][0]
                else:
                    product = ""
                gene_id = "GENE_" + feature.qualifiers["ID"][0]
                gene_dict[gene_id] = (gene, product)
    return gene_dict


gene_dict = parse_gff()
dfs = []
for d in listdir("/work/FAC/FBM/DMF/smitri/default/leakage/data/isolate_sequencing"):
    f = join(
        "/work/FAC/FBM/DMF/smitri/default/leakage/data/isolate_sequencing",
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
out.insert(12, "product", None)
out.index = range(len(out))
for i in out.index:
    gene_id = out.loc[i, "gene"]
    # If statement for annotation scenarios where non_coding_transcript_variant is annoted as gene name
    if gene_id is np.nan:
        gene, product = "", ""
    elif gene_id[:4] == "GENE":
        gene = gene_dict[gene_id][0]
        product = gene_dict[gene_id][1]
    else:
        gene = gene_id
        product = gene_id
    out.at[i, "gene"] = gene
    out.at[i, "product"] = product

for sample in set(out["sample"]):
    df = out[out["sample"] == sample]
    filter = (df["freq"] >= 0.95) & (df["depth"] >= 100)
    df = df[filter]
    df.to_csv(join("filtered_variants", sample + ".filtered.csv"), index=False)

for sample in set(out["sample"]):
    df = out[out["sample"] == sample]
    df.to_csv(join("all_variants", sample + ".csv"), index=False)


filter = (out["freq"] >= 0.95) & (out["depth"] >= 100)
out[filter].to_csv(join("filtered_variants", "all_samples.filtered.csv"), index=False)

out.to_csv(join("all_variants", "all_samples.csv"), index=False)
