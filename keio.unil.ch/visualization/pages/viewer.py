import dash
from dash import dcc, html, callback, Output, Input
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import numpy as np
from scipy.interpolate import make_smoothing_spline
import pandas as pd
import bz2
import _pickle as cPickle

def load_compressed_pickle(filename):
    data = bz2.BZ2File(filename + '.pbz2', 'rb')
    data = cPickle.load(data)
    return data

data_dict = load_compressed_pickle('data/data_dict_keio_all_strains_with_std')

all_strains = load_compressed_pickle("data/keio_non_wt_strains_list")
all_metabolites = load_compressed_pickle("data/keio_all_metabolites_list")

## load mz_
ionMz_annotation_fn = 'data/H_ionMz_annotation.csv'
df_ionMz = pd.read_csv(ionMz_annotation_fn)
df_ionMz.drop_duplicates(subset=['ionMz'], inplace=True)

mutant_color = "#0D5379"
wt_color = "#991717"

distance_fn = 'data/K_timecurve_comparison_TIC_norm_lam_10_AUC OD.csv'
dis_df = pd.read_csv(distance_fn, index_col=0)
dis_mz = dis_df.merge(df_ionMz, on='ionMz', how='left')

max_pval = 0.005
min_effect = 3
sel_idx = (dis_mz['Rel. distance'].abs() > min_effect)&(dis_mz['p-value rel'] < max_pval)
sel_df = dis_mz.loc[sel_idx].copy()

default_bar_figure = go.Figure(layout=go.Layout(
        margin=dict(l=0, r=50, t=50, b=10),
        xaxis=dict(title="Effect"),
    ))
default_bar_figure.update_xaxes(range=(-15,15))
default_bar_figure.update_yaxes(showticklabels=False)


def plot_bar_1(selected_strain,fig_layout,selected_metabolite=None,showMetaboliteFlag=False):
    idx = sel_df.Strain == selected_strain
    plot_data = sel_df.loc[idx].copy()
    # print(plot_data.head())
    plot_data.sort_values(by='Rel. distance', ascending=False, inplace=True)

    xdata = plot_data['Rel. distance'][::-1]
    ydata = plot_data['Metabolite'][::-1]

    
    fig = go.Figure(layout=fig_layout)
    fig.add_trace(
        go.Bar(
            y=ydata,
            x=xdata,
            marker_color="#0D5379",
            name='Effect size',
            width=0.9,
            orientation='h',
        )
    )
    if showMetaboliteFlag:
        for y_cat, x_val, label in zip(ydata, xdata, ydata):
            if x_val >= 0:
                # Right side of positive bar
                fig.add_annotation(
                    x=-0.5,
                    y=y_cat,
                    text=label,
                    showarrow=False,
                    xanchor="right",
                    yanchor="middle",
                    xshift=5  # pixels to the right
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
                    xshift=-5  # pixels to the left
                )
    
    fig.add_vline(
        x=0,
        line_width=1,
        line_dash="dash",
        line_color="black",
    )
    fig.update_layout(
        yaxis_title="Metabolite",
        xaxis_title="Effect",
    )
    if(len(plot_data['Metabolite'])>25):
        fig.update_yaxes(range=[len(plot_data['Metabolite'])-0.5, -0.5])
    else:
        fig.update_yaxes(range=[25-0.5, -0.5])
    min_y,max_y = np.min(plot_data['Rel. distance']), np.max(plot_data['Rel. distance'])
    range_min = min(-15,min_y)
    range_max = max(15,max_y)
    fig.update_xaxes(range=(range_min,range_max))

    fig.update_yaxes(showticklabels=False)

    if (selected_metabolite is not None) and (selected_metabolite != "Select Metabolite"):
        fig.add_shape(
            type="rect",
            x0 = range_min,
            y0 = np.where(ydata.values == selected_metabolite)[0][0]-0.5,
            x1 = range_max,
            y1 = np.where(ydata.values == selected_metabolite)[0][0]+0.5,
            line=dict(color=wt_color,width=2,dash="dashdot"),
        )


    return fig

def plot_bar_2(selected_metabolite,fig_layout,selected_strain=None):
    idx = sel_df.Metabolite == selected_metabolite
    plot_data = sel_df.loc[idx].copy()
    plot_data.sort_values(by='Rel. distance', ascending=False, inplace=True)

    xdata = plot_data['Rel. distance'][::-1]
    ydata = plot_data['Strain'][::-1]
    
    fig = go.Figure(layout=fig_layout)
    fig.add_trace(
        go.Bar(
            y=ydata,
            x=xdata,
            marker_color="#0D5379",
            name='Effect size',
            width=0.9,
            orientation='h',
        )
    )
    fig.add_vline(
        x=0,
        line_width=1,
        line_dash="dash",
        line_color="black",
    )
    fig.update_layout(
        yaxis_title="Strain",
        xaxis_title="Effect",
    )
    if(len(plot_data['Strain'])>25):
        fig.update_yaxes(range=[len(plot_data['Strain'])-0.5, -0.5])
    else:
        fig.update_yaxes(range=[25-0.5, -0.5])

    for y_cat, x_val, label in zip(ydata, xdata, ydata):
        if x_val >= 0:
            # Right side of positive bar
            fig.add_annotation(
                x=-0.5,
                y=y_cat,
                text=label,
                showarrow=False,
                xanchor="right",
                yanchor="middle",
                xshift=5  # pixels to the right
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
                xshift=-5  # pixels to the left
            )

    min_y,max_y = np.min(plot_data['Rel. distance']), np.max(plot_data['Rel. distance'])
    range_min = min(-15,min_y)
    range_max = max(15,max_y)
    fig.update_xaxes(range=(range_min,range_max))

    fig.update_yaxes(showticklabels=False)

    if(selected_strain is not None) and (selected_strain != "Select Strain"):
        fig.add_shape(
            type="rect",
            x0 = range_min,
            y0 = np.where(ydata.values == selected_strain)[0][0]-0.5,
            x1 = range_max,
            y1 = np.where(ydata.values == selected_strain)[0][0]+0.5,
            line=dict(color=wt_color,width=2,dash="dashdot"),
        )

    return fig

def plot_function(selected_strain,selected_metabolite,fig_layout):
    combined_data = data_dict[selected_strain][selected_metabolite]    
    if(combined_data is None):
        fig = go.Figure(layout=fig_layout)
        fig.add_annotation(
            x=0.5,
            y=0.5,
            text="No data available for the selected strain and metabolite.",
            showarrow=False,
            font=dict(size=16),
        )
        return fig
    
    mutant_data, wt_data = combined_data
    t_common =  np.linspace(max(mutant_data[0].min(), wt_data[0].min()), min(mutant_data[0].max(), wt_data[0].max()), 1000, endpoint=True)
    mutant_spline = make_smoothing_spline(mutant_data[0], mutant_data[1], w = 1/np.square(mutant_data[2]), lam = 10)
    wt_spline = make_smoothing_spline(wt_data[0], wt_data[1], w=1/np.square(wt_data[2]), lam = 10)

    

    fig = go.Figure(layout=fig_layout)
    fig.add_trace(
        go.Scatter(
            x=mutant_data[0],
            y=mutant_data[1],
            mode='markers',
            marker=dict(color=mutant_color,size=12,opacity=0.6,line=dict(width=1, color='DarkSlateGrey')),
            name = "Mutant Strain",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=wt_data[0],
            y=wt_data[1],
            mode='markers',
            marker=dict(color=wt_color,size=12,opacity=0.6,line=dict(width=1, color='DarkSlateGrey')),
            name = "Wild Type"
        )
    )
    fig.add_trace(
        go.Scatter(
            x=t_common,
            y=mutant_spline(t_common),
            mode='lines',
            line=dict(color=mutant_color),
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=t_common,
            y=wt_spline(t_common),
            mode='lines',
            line=dict(color=wt_color),
            showlegend=False,
        )
    )
    fig.update_layout(legend=dict(
    yanchor="top",
    y=0.99,
    xanchor="left",
    x=0.01
    ))
    return fig


# Dash layout
dash.register_page(__name__, path="/", name="Browse the data", order=1)  # '/' is home page

layout = html.Div(
    [dbc.Row([   
        # Column for strain
        dbc.Col(
            [   
                dbc.Label("Strain Selection"),
                dcc.Dropdown(
                        options=all_strains,
                        placeholder="Select strain",
                        id="strain-dropdown",
                        multi=False,
                    ),
                dbc.Checklist(
                    options={"label":"Show metabolite names"},
                    value=["label"],
                    id="metabolite-name-flag"
                ),
                dbc.Row(
                    dcc.Graph(
                        figure=default_bar_figure,
                        id="bar-1",
                    ),
                    style={"height":"80vh"},
                )           
            ],
            width=4,
        ),

        # Column for metabolite
        dbc.Col(
            [   
                dbc.Label("Metabolite Selection"),
                dcc.Dropdown(
                        options=all_metabolites,
                        placeholder="Select metabolite",
                        id="metabolite-dropdown",
                        multi=False,
                    ),
                dbc.Row(                    
                        dcc.Graph(
                            figure=default_bar_figure,
                            id="bar-2",
                        ),
                        style={"height":"80vh"},                    
                ),                    
            ],
            width=4,
        ),

        # Column for plot
         dbc.Col(
                dcc.Graph(
                    figure={},
                    id="controls-and-graph",
                    style={"width":"100%", "height":"40vh"}
                ),
                width=4,
                style={"height":"80vh","justify-content": "center", "align-items": "center", "display": "flex","padding-left":"2vw"},
            ),
        

    ],),
    dbc.Col(
            html.Details(
                [
                    html.Summary(html.Span("Metadata")),
                    html.Span(
                        children=[""],
                        id="metadata-table",
                    ),
                ],
            ),
            width=11,
    ),],
)


# Callback for plotting the selected data
@callback(
    [
        Output(component_id="controls-and-graph", component_property="figure"),
        Output(component_id="metadata-table", component_property="children"),
    ],
    [
        Input("strain-dropdown", "value"),
        Input("metabolite-dropdown", "value"),
    ],
)

        # Output(component_id="bar-2", component_property="figure"),
def update_graph_view(
    chosen_strains, chosen_metabolites
):
    fig_layout = go.Layout(
        margin=dict(l=0, r=50, t=50, b=10),
        xaxis=dict(title="AUC OD"),
        yaxis=dict(title="Median Z-score"),
    )
    
    if (
        chosen_strains == "Select Strain"
        or chosen_strains == None
        or chosen_metabolites == "Select Metabolite"
        or chosen_metabolites == None
    ):
        fig = go.Figure(layout=fig_layout)
        return fig, []

    fig = plot_function(chosen_strains, chosen_metabolites, fig_layout)

    return (
        fig, []
    )

# Callback for plotting the selected data
@callback(
        Output(component_id="bar-1", component_property="figure"),
        Input("strain-dropdown", "value"),
        Input("metabolite-dropdown", "value"),
        Input("metabolite-name-flag", "value"),
)
        
        
def update_bar_graph1(chosen_strains, chosen_metabolites, show_metabolite_flag):
    fig_layout = go.Layout(
        margin=dict(l=0, r=50, t=50, b=10),
        xaxis=dict(title="Effect"),
    )
    if (
        chosen_strains == "Select Strain"
        or chosen_strains == None
    ):
        return default_bar_figure
    fig = plot_bar_1(chosen_strains,fig_layout=fig_layout,selected_metabolite=chosen_metabolites,showMetaboliteFlag=show_metabolite_flag)

    return fig

# Callback for plotting the selected data
@callback(
        Output(component_id="bar-2", component_property="figure"),
        Input("metabolite-dropdown", "value"),
        Input("strain-dropdown", "value"),
)

def update_bar_graph2(chosen_metabolites, chosen_strains):
    fig_layout = go.Layout(
        margin=dict(l=0, r=50, t=50, b=10),
        xaxis=dict(title="Effect"),
    )
    if (
        chosen_metabolites == "Select Metabolite"
        or chosen_metabolites == None
    ):
        return default_bar_figure
    fig = plot_bar_2(chosen_metabolites,fig_layout=fig_layout,selected_strain=chosen_strains)

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