#!/usr/bin/env python
#
# Reserve 1 CPUs for this job
#
# SBATCH --cpus-per-task=1
# SBATCH --mem=2G
#
# Request it to run this for DD:HH:MM with ?G per core
#
# SBATCH --time=24:00:00
#
import subprocess
from os.path import join
import os

data = "/work/FAC/FBM/DMF/smitri/evomicrocomm/seq_snorre/data/meta_sequencing"


def submit(files):
    """Basic snakemake calling taking the desired output file as input."""
    cluster_config = '--cluster-config cluster.json --cluster \
        "sbatch --mem={cluster.mem} -t {cluster.time} -c {threads}"'
    cmd = [
        "snakemake",
        "--latency-wait",
        "5",
        "--rerun-incomplete",
        "-j",
        "500",
        cluster_config,
        files,
    ]
    subprocess.call(" ".join(cmd), shell=True)


def sample_caller(output_file):
    output = [sample for sample in os.listdir(data)]
    files = join(data, "{" + ",".join(output) + "}", output_file)
    submit(files)


sample_caller("var.annotated.vcf")
