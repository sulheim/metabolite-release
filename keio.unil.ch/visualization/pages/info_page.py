import dash
from dash import dcc, html, callback
import markdown

# Load markdown content
def load_markdown_content(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    return content

# Register the new page
dash.register_page(__name__, path="/info", name="Information", order=2)

layout = html.Div(
    [
        # html.H2("Information Page"),
        html.Img(src='data/illustration.png', style={'width': '100%', 'height': 'auto'}),
        dcc.Markdown(load_markdown_content('data/info.md'))
    ]
)