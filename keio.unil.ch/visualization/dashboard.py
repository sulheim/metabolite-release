from os.path import join

import dash
from dash import html, dcc, Input, Output, callback, ctx
import dash_bootstrap_components as dbc

app = dash.Dash(
    __name__,
    use_pages=True,
    external_stylesheets=[dbc.themes.BOOTSTRAP, join("assets", "style_new.css")],
    suppress_callback_exceptions=True,
)
server = app.server

pages = list(dash.page_registry.values())
pages = sorted(pages, key=lambda p: p.get("order", 0))
default_path = pages[0]["path"] if pages else "/"
valid_paths = {p["path"] for p in pages}

tabs = dcc.Tabs(
    id="page-tabs",
    value=default_path,
    children=[dcc.Tab(label=p["name"], value=p["path"]) for p in pages],
    persistence=True,
    persistence_type="session",
)

app.layout = dbc.Container(
    [
        dcc.Location(id="url"),
        dbc.Row(
            html.H1(
                "KEIO Metabolite Release Dashboard",
                className="pt-4 ps-5",
            )
        ),
        dbc.Row(
            tabs,
            class_name="ps-5 pt-3",
        ),
        dbc.Row(
            dash.page_container,
            class_name="ps-5 pt-4",
        ),
    ],
    fluid=True,
)


# Callbacks to call pages if tab is clicked
@callback(
    Output("url", "pathname"),
    Output("page-tabs", "value"),
    Input("url", "pathname"),
    Input("page-tabs", "value"),
)
def sync_tabs_and_url(pathname, tab_value):
    trigger = ctx.triggered_id

    if trigger == "page-tabs":
        path = tab_value if tab_value in valid_paths else default_path
        return path, path

    path = pathname if pathname in valid_paths else default_path
    return path, path


if __name__ == "__main__":
    app.run(debug=True)
