import dash
from dash import dcc, html, callback


# Load markdown content
def load_markdown_content(filepath):
    with open(filepath, "r") as f:
        content = f.read()
    return content


# Register the new page
dash.register_page(__name__, path="/info", name="About", order=2)

layout = html.Div(
    [
        # html.H2("Information Page"),
        html.Div(
            html.Img(
                src="assets/illustration.png", style={"width": "90%", "maxWidth": "600px", "height": "auto"}
            ),
            style={"textAlign": "center"},
        ),
        dcc.Markdown(load_markdown_content("data/info.md")),
    ],
    style={
        "justify-content": "center",
        "align-items": "center",
        "display": "flex",
        "flex-direction": "column",
    },
)
