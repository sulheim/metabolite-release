import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
from style import *
from plotly.subplots import make_subplots
from os import listdir
from os.path import join


dir = "/work/FAC/FBM/DMF/smitri/evomicrocomm/seq_snorre/data/meta_sequencing"
dirs = sorted([d for d in listdir(dir)])


def fig1_coverage():
    row, col = 1, 1
    fig = make_subplots(rows=2, cols=3, subplot_titles=dirs)
    for i, d in enumerate(dirs):
        print(d)
        f = join(dir, d, "coverage.txt")
        df = pd.read_csv(f, sep="\t", names=["chrom", "pos", "depth"])
        fig.add_trace(
            go.Scatter(
                x=df["pos"][::100],
                y=df["depth"][::100],
                mode="markers",
                marker=dict(color=colors["blue"]),
                showlegend=False,
            ),
            row=row,
            col=col,
        )
        fig.update_xaxes(title="Position"), fig.update_yaxes(
            title="Coverage", type="log"
        )

        col += 1
        if i == 2:
            row = 2
            col = 1

    fig = style_plot(fig, marker_size=1)
    fig.write_image("plots/fig1_coverage.svg")


def fig2_allele_frequency():
    df = pd.read_csv("results/freebayes_snps/all_samples.filtered.csv")
    df = df.sort_values(by="sample")
    color_map = {dir: colors["blue"] for dir in dirs}
    fig = px.strip(x=df["sample"], y=df["freq"], color_discrete_map=color_map)
    fig.update_xaxes(title="Sample"), fig.update_yaxes(title="Allele frequency")
    fig = style_plot(fig, marker_size=3)
    fig.write_image("plots/fig2_frequencies.svg")


fig2_allele_frequency()


def fig3_total_allele_frequency():
    fig = go.Figure()
    for d in dirs:
        f = "results/freebayes_snps/" + d + ".filtered.csv"
        df = pd.read_csv(f)
        fig.add_trace(
            go.Scatter(
                x=[d],
                y=[sum(df["freq"])],
                showlegend=False,
                marker=dict(color=colors["blue"]),
            )
        )
    fig = style_plot(fig, marker_size=5)
    fig.update_layout(width=width, height=height)
    fig.write_image("plots/fig3_total_freq.svg")
