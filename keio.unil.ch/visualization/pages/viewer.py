import dash
from dash import dcc, html, callback, Output, Input, dash_table
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
from scipy.interpolate import make_smoothing_spline
import pandas as pd
import bz2
import _pickle as cPickle

# Colors
colormap = px.colors.qualitative.Set2
mutant_color = colormap[1]  # "#0D5379"
wt_color = colormap[0]  # "#991717"
selected_color = colormap[1]
bar_color = colormap[2]

# Label size
label_dict = dict(
    # family="Arial, sans-serif",
    size=14,  # Global font size
    color="black",
)

default_font_size = 14

# Starting point for data loading
first_strain = "putP"
first_metabolite = "Proline"
strain_placeholder = "-"
metabolite_placeholder = "-"


def load_compressed_pickle(filename):
    data = bz2.BZ2File(filename + ".pbz2", "rb")
    data = cPickle.load(data)
    return data


data_dict = load_compressed_pickle("data/new_data_dict_keio_all_strains_with_std")

all_strains = load_compressed_pickle("data/keio_non_wt_strains_list")
all_metabolites = load_compressed_pickle("data/keio_all_metabolites_list")

# Add delta
strain_options = [{"label": f"Δ{strain}", "value": strain} for strain in all_strains]

## load metabolite metadata
ionMz_annotation_fn = "data/H_ionMz_annotation.csv"
df_ionMz = pd.read_csv(ionMz_annotation_fn)
df_ionMz.drop_duplicates(subset=["ionMz"], inplace=True)

# Load estimated effects
distance_fn = "data/K_timecurve_comparison_TIC_norm_lam_10_AUC OD.csv"
dis_df = pd.read_csv(distance_fn, index_col=0)
dis_mz = dis_df.merge(df_ionMz, on="ionMz", how="left")

# Load sample metadata
# sample_metadata_fn = 'data/I_sample_metadata_keio.csv'
# sample_metadata_df = pd.read_csv(sample_metadata_fn, index_col=0)

# Load transporter metadata
transporter_metadata_fn = "data/C_transporters_short_metadata.csv"
transporter_metadata_df = pd.read_csv(transporter_metadata_fn)


max_pval = 0.005
min_effect = 3
sel_idx = (dis_mz["Rel. distance"].abs() > min_effect) & (
    dis_mz["p-value rel"] < max_pval
)
sel_df = dis_mz.loc[sel_idx].copy()

default_bar_figure = go.Figure(
    layout=go.Layout(
        margin=dict(l=0, r=50, t=50, b=10),
        xaxis=dict(title="Effect"),
    )
)
default_bar_figure.update_xaxes(range=(-15, 15))
default_bar_figure.update_yaxes(showticklabels=False)


def plot_bar_1(
    selected_strain, fig_layout, selected_metabolite=None, showMetaboliteFlag=False
):
    idx = sel_df.Strain == selected_strain
    plot_data = sel_df.loc[idx].copy()
    # print(plot_data.head())
    plot_data.sort_values(by="Rel. distance", ascending=False, inplace=True)

    xdata = plot_data["Rel. distance"][::-1]
    ydata = plot_data["Metabolite"][::-1]

    # Default color for all bars
    colors = [bar_color] * len(ydata)
    # Highlight selected metabolite
    if (selected_metabolite is not None) and (selected_metabolite != "-"):
        selected_idx = np.where(ydata.values == selected_metabolite)[0]
        if len(selected_idx):
            colors[selected_idx[0]] = (
                selected_color  # "#991717"  # Use wt_color or any highlight color
            )

    fig = go.Figure(layout=fig_layout)
    # fig.
    fig.add_trace(
        go.Bar(
            y=ydata,
            x=xdata,
            marker_color=colors,
            name="Effect size",
            width=0.9,
            orientation="h",
        )
    )
    # Truncate to 20 characters
    max_length = 20
    truncated_labels = [
        label[:max_length] + "..." if len(label) > max_length + 3 else label
        for label in ydata
    ]

    if showMetaboliteFlag:
        for y_cat, x_val, label in zip(ydata, xdata, truncated_labels):
            if x_val >= 0:
                # Right side of positive bar
                fig.add_annotation(
                    x=-0.5,
                    y=y_cat,
                    text=label,
                    showarrow=False,
                    xanchor="right",
                    yanchor="middle",
                    xshift=5,  # pixels to the right
                )
            else:
                # Left side of negative bar
                fig.add_annotation(
                    x=0.5,
                    y=y_cat,
                    text=label,
                    showarrow=False,
                    xanchor="left",
                    yanchor="middle",
                    xshift=-5,  # pixels to the left
                )

    fig.add_vline(
        x=0,
        line_width=1,
        # line_dash="dash",
        line_color="black",
    )
    fig.add_annotation(
        x=0.01,
        y=1,
        text="Click on a bar to select that metabolite",
        xref="paper",
        yref="paper",
        xanchor="left",
        yanchor="bottom",
        # bgcolor = "white",
        showarrow=False,
        font=dict(size=default_font_size, color="gray"),
    )
    fig.update_layout(yaxis_title="Metabolite", xaxis_title="Effect", font=label_dict)

    if len(plot_data["Metabolite"]) > 25:
        fig.update_yaxes(range=[len(plot_data["Metabolite"]) - 0.5, -0.5])
    else:
        fig.update_yaxes(range=[25 - 0.5, -0.5])
    min_y, max_y = np.min(plot_data["Rel. distance"]), np.max(
        plot_data["Rel. distance"]
    )
    range_min = min(-15, min_y)
    range_max = max(15, max_y)
    fig.update_xaxes(range=(range_min, range_max))
    fig.update_yaxes(showticklabels=False)

    return fig


def plot_bar_2(selected_metabolite, fig_layout, selected_strain=None):
    idx = sel_df.Metabolite == selected_metabolite
    plot_data = sel_df.loc[idx].copy()
    plot_data.sort_values(by="Rel. distance", ascending=False, inplace=True)

    xdata = plot_data["Rel. distance"][::-1]
    ydata = plot_data["Strain"][::-1]

    # Default color for all bars
    colors = [bar_color] * len(ydata)
    # Highlight selected metabolite
    if (selected_strain is not None) and (selected_strain != "-"):
        selected_idx = np.where(ydata.values == selected_strain)[0]
        if len(selected_idx):
            colors[selected_idx[0]] = (
                selected_color  # "#991717"  # Use wt_color or any highlight color
            )

    fig = go.Figure(layout=fig_layout)
    fig.add_trace(
        go.Bar(
            y=ydata,
            x=xdata,
            marker_color=colors,
            name="Effect size",
            width=0.9,
            orientation="h",
        )
    )
    fig.add_vline(
        x=0,
        line_width=1,
        # line_dash="dash",
        line_color="black",
    )
    fig.update_layout(yaxis_title="Strain", xaxis_title="Effect", font=label_dict)
    if len(plot_data["Strain"]) > 25:
        fig.update_yaxes(range=[len(plot_data["Strain"]) - 0.5, -0.5])
    else:
        fig.update_yaxes(range=[25 - 0.5, -0.5])

    ytext = [f"Δ{strain}" for strain in ydata]
    for y_cat, x_val, label in zip(ydata, xdata, ytext):
        if x_val >= 0:
            # Right side of positive bar
            fig.add_annotation(
                x=-0.5,
                y=y_cat,
                text=label,
                showarrow=False,
                xanchor="right",
                yanchor="middle",
                xshift=5,  # pixels to the right
            )
        else:
            # Left side of negative bar
            fig.add_annotation(
                x=0.5,
                y=y_cat,
                text=label,
                showarrow=False,
                xanchor="left",
                yanchor="middle",
                xshift=-5,  # pixels to the left
            )
    fig.add_annotation(
        x=0.01,
        y=1,
        text="Click on a bar to select that KO mutant",
        xref="paper",
        yref="paper",
        xanchor="left",
        yanchor="bottom",
        # bgcolor = "white",
        showarrow=False,
        font=dict(size=default_font_size, color="gray"),
    )
    min_y, max_y = np.min(plot_data["Rel. distance"]), np.max(
        plot_data["Rel. distance"]
    )
    range_min = min(-15, min_y)
    range_max = max(15, max_y)
    fig.update_xaxes(range=(range_min, range_max))

    fig.update_yaxes(showticklabels=False)

    return fig


def plot_function(
    selected_strain, selected_metabolite, fig_layout, xaxis_type="AUC OD"
):
    combined_data = data_dict[selected_strain][selected_metabolite]
    if combined_data is None:
        fig = go.Figure(layout=fig_layout)
        fig.add_annotation(
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            text="No data available for the selected strain and metabolite.",
            showarrow=False,
            font=dict(size=16),
        )
        return fig

    mutant_data, wt_data = combined_data
    zi = 1
    zsi = 3
    if xaxis_type == "Time [h]":
        xi = 0
    else:
        xi = 2

    wt_data = wt_data[:, np.argsort(wt_data[xi])]
    mutant_data = mutant_data[:, np.argsort(mutant_data[xi])]

    t_m = np.linspace(mutant_data[xi].min(), mutant_data[xi].max(), 100, endpoint=True)
    t_wt = np.linspace(wt_data[xi].min(), wt_data[xi].max(), 100, endpoint=True)

    mutant_spline = make_smoothing_spline(
        mutant_data[xi], mutant_data[zi], w=1 / np.square(mutant_data[zsi]), lam=10
    )
    wt_spline = make_smoothing_spline(
        wt_data[xi], wt_data[zi], w=1 / np.square(wt_data[zsi]), lam=10
    )

    fig = go.Figure(layout=fig_layout)
    fig.add_trace(
        go.Scatter(
            x=mutant_data[xi],
            y=mutant_data[zi],
            mode="markers",
            marker=dict(
                color=mutant_color,
                size=12,
                opacity=0.6,
                line=dict(width=1, color="DarkSlateGrey"),
            ),
            name=f"Δ{selected_strain}",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=wt_data[xi],
            y=wt_data[zi],
            mode="markers",
            marker=dict(
                color=wt_color,
                size=12,
                opacity=0.6,
                line=dict(width=1, color="DarkSlateGrey"),
            ),
            name="Wild Type",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=t_m,
            y=mutant_spline(t_m),
            mode="lines",
            line=dict(color=mutant_color),
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=t_wt,
            y=wt_spline(t_wt),
            mode="lines",
            line=dict(color=wt_color),
            showlegend=False,
        )
    )
    fig.update_layout(
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
            font={"size": default_font_size},
        ),
        font=label_dict,
    )
    return fig


def metadata_function(selected_strain, selected_metabolite):
    # Get transporter info
    mutant_info = transporter_metadata_df.loc[
        transporter_metadata_df["Gene Name"] == selected_strain
    ]
    if mutant_info.empty:
        return html.Span("No metadata available.")
    mutant_info = mutant_info.iloc[0]

    # Get metabolite info
    metabolite_info = df_ionMz.loc[df_ionMz["Metabolite"] == selected_metabolite]
    if metabolite_info.empty:
        return html.Span("No metadata available.")
    metabolite_info = metabolite_info.iloc[0]

    # Transporter metadata table
    transporter_df = pd.DataFrame(
        {
            "Property": [
                "Gene",
                "JW ID",
                "Annotation",
                "Location",
                "Type",
            ],
            "Value": [
                selected_strain,
                mutant_info["JW ID"],
                mutant_info["Annotation"].capitalize(),
                mutant_info["Location"].capitalize(),
                mutant_info["Type"],
            ],
        }
    )

    # Metabolite metadata table
    metabolite_df = pd.DataFrame(
        {
            "Property": [
                "Name",
                "Ion m/z",
                "All KEGG IDs",
                "Other annotations",
                "Class",
            ],
            "Value": [
                selected_metabolite,
                metabolite_info["ionMz"],
                metabolite_info["KEGG"],
                metabolite_info["Alt. annotations"],
                metabolite_info["Classification"],
            ],
        }
    )

    return dbc.Row(
        [
            dbc.Label(
                "Metadata",
                style={
                    "fontWeight": "regular",
                    "fontSize": "18px",
                    "marginBottom": "0vh",
                },
            ),
            dbc.Col(
                dash_table.DataTable(
                    data=transporter_df.to_dict("records"),
                    columns=[
                        {"name": "KO mutant", "id": "Property"},
                        {"name": "", "id": "Value"},
                    ],
                    style_cell={"textAlign": "left", "padding": "2px"},
                    style_header={
                        "backgroundColor": mutant_color,
                        "fontWeight": "bold",
                    },
                    style_as_list_view=True,
                    style_data={"fontSize": f"{default_font_size}px"},
                    style_cell_conditional=[
                        {"if": {"column_id": "Property"}, "width": "30%"},
                        # {'if': {'column_id': 'Value'}, 'width': '200px'},
                    ],
                ),
                width=6,
            ),
            dbc.Col(
                dash_table.DataTable(
                    data=metabolite_df.to_dict("records"),
                    columns=[
                        {"name": "Metabolite", "id": "Property"},
                        {"name": "", "id": "Value"},
                    ],
                    style_cell={"textAlign": "left", "padding": "2px"},
                    style_header={
                        "backgroundColor": bar_color,
                        "fontWeight": "bold",
                        "fontSize": f"{default_font_size}px",
                    },
                    style_data={"fontSize": f"{default_font_size}px"},
                    style_cell_conditional=[
                        {"if": {"column_id": "Property"}, "width": "40%"},
                        # {'if': {'column_id': 'Value'}, 'width': '200px'},
                    ],
                    style_as_list_view=True,
                ),
                width=6,
            ),
        ],
        style={"marginTop": "1vh", "marginBottom": "1vh", "marginRight": "2vw"},
    )


def plot_OD(selected_strain, fig_layout):

    if (selected_strain is None) or (selected_strain == "-"):
        fig = go.Figure(layout=fig_layout)
        fig.add_annotation(
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            text="No KO mutant selected.",
            showarrow=False,
            font=dict(size=default_font_size),
        )
        return fig

    mutant_data, wt_data = data_dict[selected_strain]["OD"]
    xi = 0
    yi = 1
    ysi = 2

    fig = go.Figure(layout=fig_layout)
    fig.add_trace(
        go.Scatter(
            x=mutant_data[xi],
            y=mutant_data[yi],
            mode="markers+lines",
            error_y=dict(type="data", array=mutant_data[ysi], visible=True),
            marker=dict(
                color=mutant_color,
                size=12,
                opacity=0.6,
                line=dict(width=1, color="DarkSlateGrey"),
            ),
            name=f"Δ{selected_strain}",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=wt_data[xi],
            y=wt_data[yi],
            mode="markers+lines",
            error_y=dict(type="data", array=wt_data[ysi], visible=True),
            marker=dict(
                color=wt_color,
                size=12,
                opacity=0.6,
                line=dict(width=1, color="DarkSlateGrey"),
            ),
            name="Wild Type",
        )
    )

    fig.update_layout(
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
            font={"size": default_font_size},
        ),
        font=label_dict,
    )
    return fig


# Dash layout
dash.register_page(__name__, path="/", name="Browse data", order=1)  # '/' is home page

layout = html.Div(
    [
        dbc.Row(
            [
                # Column for strain
                dbc.Col(
                    [
                        dbc.Label(
                            "Select transporter KO mutant",
                        ),
                        dcc.Dropdown(
                            options=strain_options,
                            placeholder="-",
                            id="strain-dropdown",
                            value=first_strain,
                            multi=False,
                            className="dropdownS",
                            style={"marginBottom": "0vh"},
                        ),
                        dbc.Checklist(
                            options={"label": "Show metabolite names"},
                            value=["label"],
                            id="metabolite-name-flag",
                            style={
                                "marginTop": "0vh",
                                "marginBottom": "0vh",
                                "height": "5vh",
                            },
                        ),
                        # dbc.Row(
                        dcc.Graph(
                            figure=default_bar_figure,
                            id="bar-1",
                            style={"height": "80vh", "marginTop": "0vh"},
                        ),
                        # )
                    ],
                    width=3,
                ),
                # Column for metabolite
                dbc.Col(
                    [
                        dbc.Label("Metabolite Selection"),
                        dcc.Dropdown(
                            options=all_metabolites,
                            placeholder="-",
                            id="metabolite-dropdown",
                            value=first_metabolite,
                            multi=False,
                            className="dropdownM",
                        ),
                        dbc.Row(
                            dcc.Graph(
                                figure=default_bar_figure,
                                id="bar-2",
                            ),
                            style={"height": "80vh", "marginTop": "5vh"},
                        ),
                    ],
                    width=3,
                ),
                # Column for plot and metadata
                dbc.Col(
                    [
                        dbc.Label("Select x-axis"),
                        dcc.Dropdown(
                            options=["Time [h]", "AUC OD"],
                            value="AUC OD",  # default value
                            id="xaxis-dropdown",
                            className="dropdownM",
                            clearable=False,
                        ),
                        dcc.Graph(
                            figure={},
                            id="controls-and-graph",
                            style={
                                "width": "100%",
                                "height": "30vh",
                                "marginTop": "2.5vh",
                            },
                        ),
                        dcc.Graph(
                            figure={},
                            id="od-graph",
                            style={
                                "width": "100%",
                                "height": "30vh",
                                "marginTop": "0vh",
                            },
                        ),
                        # html.Details(
                        #     [
                        # html.Summary(html.Span("Metadata")),
                        html.Span(
                            children=[""],
                            id="metadata-table",
                            style={"marginTop": "1vh", "marginBottom": "1vh"},
                        ),
                        # ],
                        # style={"marginTop": "0vh", "marginBottom": '1vh'}  #
                        # ),
                    ],
                    width=6,
                    style={"height": "80vh", "padding-left": "2vw"},
                ),
            ],
            style={"marginLeft": "1vw"},
        ),
        # Remove the separate dbc.Col for metadata here
    ]
)


# Callback for plotting the selected data
@callback(
    Output(component_id="controls-and-graph", component_property="figure"),
    # Output(component_id="metadata-table", component_property="children"),
    [
        Input("strain-dropdown", "value"),
        Input("metabolite-dropdown", "value"),
        Input("xaxis-dropdown", "value"),
    ],
)


# Output(component_id="bar-2", component_property="figure"),
def update_graph_view(chosen_strains, chosen_metabolites, xaxis_type):
    fig_layout = go.Layout(
        margin=dict(l=0, r=50, t=50, b=10),
        xaxis=dict(title=xaxis_type),
        yaxis=dict(title="Median Z-score"),
    )

    if (
        chosen_strains == "Select Strain"
        or chosen_strains == None
        or chosen_metabolites == "-"
        or chosen_metabolites == None
    ):
        fig = go.Figure(layout=fig_layout)
        return fig

    fig = plot_function(
        chosen_strains, chosen_metabolites, fig_layout, xaxis_type=xaxis_type
    )
    return fig


# Plot OD
@callback(
    Output(component_id="od-graph", component_property="figure"),
    Input("strain-dropdown", "value"),
)
def update_od_graph(selected_strain):
    fig_layout = go.Layout(
        margin=dict(l=0, r=50, t=50, b=10),
        xaxis=dict(title="Time [h]"),
        yaxis=dict(title="OD600"),
    )
    if (selected_strain == None) or (selected_strain == "-"):
        fig = go.Figure(layout=fig_layout)
        return fig
    fig = plot_OD(selected_strain, fig_layout)
    return fig


# Callback for plotting the selected data
@callback(
    Output(component_id="bar-1", component_property="figure"),
    Input("strain-dropdown", "value"),
    Input("metabolite-dropdown", "value"),
    Input("metabolite-name-flag", "value"),
)
def update_bar_graph1(chosen_strains, chosen_metabolites, show_metabolite_flag):
    fig_layout = go.Layout(
        margin=dict(l=0, r=10, t=50, b=10),
        xaxis=dict(title="Effect"),
    )
    if chosen_strains == "-" or chosen_strains == None:
        return default_bar_figure
    fig = plot_bar_1(
        chosen_strains,
        fig_layout=fig_layout,
        selected_metabolite=chosen_metabolites,
        showMetaboliteFlag=show_metabolite_flag,
    )

    return fig


# Callback for plotting the selected data
@callback(
    Output(component_id="bar-2", component_property="figure"),
    Input("metabolite-dropdown", "value"),
    Input("strain-dropdown", "value"),
)
def update_bar_graph2(chosen_metabolites, chosen_strains):
    fig_layout = go.Layout(
        margin=dict(l=0, r=10, t=50, b=10),
        xaxis=dict(title="Effect"),
    )
    if chosen_metabolites == "-" or chosen_metabolites == None:
        return default_bar_figure
    fig = plot_bar_2(
        chosen_metabolites, fig_layout=fig_layout, selected_strain=chosen_strains
    )

    return fig


@callback(
    Output("metabolite-dropdown", "value"),
    Input("bar-1", "clickData"),
    prevent_initial_call=True,
)
def update_metabolite_dropdown_on_bar_click(clickData):
    if clickData is None:
        return dash.no_update
    clicked_metabolite = clickData["points"][0]["y"]
    return clicked_metabolite


@callback(
    Output("strain-dropdown", "value"),
    Input("bar-2", "clickData"),
    prevent_initial_call=True,
)
def update_strain_dropdown_on_bar_click(clickData):
    if clickData is None:
        return dash.no_update
    clicked_strain = clickData["points"][0]["y"]
    return clicked_strain


# Callback for updating metadata
@callback(
    Output("metadata-table", "children"),
    Input("strain-dropdown", "value"),
    Input("metabolite-dropdown", "value"),
)
def update_metadata_table(selected_strain, selected_metabolite):
    return metadata_function(selected_strain, selected_metabolite)
