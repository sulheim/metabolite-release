import dash
from dash import html
import dash_bootstrap_components as dbc
from os.path import join


app = dash.Dash(
    __name__,
    use_pages=True,
    external_stylesheets=[dbc.themes.BOOTSTRAP, join("assets", "style_new.css")],
)
server = app.server
sidebar = dbc.Nav(
    [
        dbc.NavLink(
            [
                html.Div(page["name"], className="ms-2"),
            ],
            href=page["path"],
            active="exact",
        )
        for page in dash.page_registry.values()
    ],
    vertical=True,
    pills=True,
    className="bg-light",
)

app.layout = dbc.Container(
    [
        dbc.Col(
            dbc.Row(
                html.H1("KEIO Metabolite Release Dashboard"),
                class_name="pt-4",
            ),
            class_name="ps-5",
        ),
        dbc.Row(
            [
                dbc.Col([sidebar], class_name="ps-5 pt-4", width=2),
                dbc.Col([dash.page_container], class_name="ps-5 pt-4", width=10),
            ]
        ),
    ],
    fluid=True,
)


if __name__ == "__main__":
    app.run(debug=True)
