width, height = 500, 400
colors = {"blue": "#000080"}


def style_plot(
    fig,
    marker_size=3,
    top_margin=20,
    left_margin=40,
    right_margin=0,
    buttom_margin=50,
    font_size=14,
    line_thickness=3,
):
    """Style function for figures setting fot size and true black color."""
    fig.update_layout(
        {
            "plot_bgcolor": "#FFFFFF",
            "paper_bgcolor": "#FFFFFF",
        },
        font={"size": font_size, "color": "black"},
    )
    for d in fig["data"]:
        d["marker"]["size"] = marker_size
        d["line"]["width"] = line_thickness
        # d["error_y"]["thickness"] = line_thickness
    for a in fig["layout"]["annotations"]:
        a["font"]["size"] = font_size
        a["font"]["color"] = "black"
    fig["layout"]["title"]["font"]["size"] = font_size
    fig["layout"]["title"]["font"]["color"] = "black"
    fig["layout"]["legend"]["title"]["font"]["size"] = font_size
    fig["layout"]["legend"]["title"]["font"]["color"] = "black"

    fig.update_layout(
        margin=dict(l=left_margin, r=right_margin, t=top_margin, b=buttom_margin),
        hoverlabel=dict(font_size=font_size),
    )
    gridline_width = 1
    fig.update_yaxes(
        title_standoff=0,
        gridcolor="black",
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=gridline_width,
    )
    fig.update_xaxes(
        title_standoff=0,
        gridcolor="black",
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=gridline_width,
    )
    fig.for_each_xaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    fig.for_each_yaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    return fig
