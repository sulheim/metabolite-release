import subprocess
from os.path import join
import os

data = "/work/FAC/FBM/DMF/smitri/default/leakage/data/isolate_sequencing"


def submit(files):
    """Basic snakemake calling taking the desired output file as input."""
    cluster_config = '--cluster-config cluster.json --cluster \
        "sbatch --mem={cluster.mem} -t {cluster.time} -c {threads}"'
    cmd = ["snakemake", "--rerun-incomplete", "-j", "500", cluster_config, files]
    subprocess.call(" ".join(cmd), shell=True)


def sample_caller(output_file):
    output = [sample for sample in os.listdir(data)]
    files = join(data, "{" + ",".join(output) + "}", output_file)
    submit(files)


sample_caller("snippy/snps.tab")
