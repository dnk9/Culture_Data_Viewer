import plotly.graph_objects as go
import pandas as pd


def plot(cdf: pd.DataFrame, x_axes_list: list) -> go.Figure:
    """
    Create a go.Figure from the input cdf, plotting the columns in the x_axes_list.

    The "domain_choice_axis_numb" is needed to specify the ratio between chart and space reserved for the y_axes names
    to be displayed. "domain" is a property of the x_axis. It accepts values between 0 and 1. The chart is plotted in
    the region defined by this two values. When multiple y_axes are needed, the space which they require should be
    kept free, by limiting the dimension of the chart.
    The "domain_choice_axis_numb" contains different values of "domain" tailored for different numbers of y_axes.


    Parameters
    ----------
    cdf: pd.DataFrame
    x_axes_list: list
        List of dicts. Each dict contains a list of columns which are plotted with the same y axis. It specifies the
        title of the axis, the font color, the tick format and the plotting mode ("lines", "lines+markers" or "markers")

    Returns
    -------

    """
    figx = go.Figure()

    domain_choice_axis_numb = {0: {"domain": [0, 1]}, 1: {"domain": [0, 1]}, 2: {"domain": [0.1, 0.0]},
                               3: {"domain": [0.10, 0.85]}}
    positions = [0, 1, 0.1, 0.9]
    sides = ["left", "right", "left", "right"]

    n_axis = len(x_axes_list)
    layout_axis = {}
    layout_axis["xaxis"] = domain_choice_axis_numb[n_axis - 1]

    for n, ax in enumerate(x_axes_list):
        yaxis = f"y{n + 1}"
        yaxis_long = f"yaxis{n + 1}"
        yaxis_title = ax["title"] if ax["title"] else ", ".join(ax["cols"])
        titlefont = ax["titlefont"]
        tickfont = ax["tickfont"]
        mode = ax["mode"]
        position = positions[n]
        side = sides[n]

        overlaying = "y" if n > 0 else None

        print(n + 1, yaxis, yaxis_title, position, side, overlaying)

        for col in ax["cols"]:
            figx.add_trace(go.Scatter(x=cdf.index, y=cdf[col], name=col, mode=mode, connectgaps=True, yaxis=yaxis))

        layout_axis[yaxis_long] = (
            dict(title=yaxis_title, position=position, overlaying=overlaying, titlefont=titlefont, tickfont=tickfont))

    figx.update_xaxes(dtick=2, rangemode="nonnegative")
    figx.update_layout(layout_axis)
    figx.update_layout(title={"text":"Culture Data Viewer",'x':0.5, "xanchor":"center"}, xaxis_title="Time (h)")

    return figx
