import pandas as pd
import plotly.express as px

stats = pd.concat(pd.read_csv(f, sep="\t") for f in snakemake.input.stats)
template = "plotly_white"

fig0 = px.scatter(
    stats,
    x="recall",
    y="precision",
    color="fdr",
    size="coverage",
    template=template,
    color_continuous_scale="viridis",
)
fig0.write_image(snakemake.output.plot0, width=800, height=400)

fig1 = px.scatter(
    stats,
    x="recall",
    y="precision",
    color="coverage",
    size="fdr",
    template=template,
    color_continuous_scale="viridis",
)
fig1.write_image(snakemake.output.plot1, width=800, height=400)

import plotly.graph_objects as go

fig = go.Figure(
    layout=go.Layout(
        template=template,
        xaxis=go.layout.XAxis(title="coverage"),
        yaxis=go.layout.YAxis(title="value"),
    )
)


def get_traces(fdr: float):
    df = stats.query(f"fdr == {fdr}")
    recall = go.Scatter(
        x=df["coverage"],
        y=df["recall"],
        name=f"recall",
        marker=dict(color="#0077bb" if fdr == 0.0 else "#33bbee"),
        mode="markers",
    )

    precision = go.Scatter(
        x=df["coverage"],
        y=df["precision"],
        name=f"precision",
        marker=dict(color="#ee7733" if fdr == 0.0 else "#cc3311"),
        mode="markers",
    )
    return precision, recall


# precision0, recall0 = get_traces(0.0)
precision1, recall1 = get_traces(1.0)
# fig.add_trace(precision0)
# fig.add_trace(recall0)
fig.add_trace(precision1)
fig.add_trace(recall1)
fig.update_layout(
    margin=dict(l=0, r=0, b=0, t=0),
    yaxis_range=[0, 1.05],
)
fig.write_image(snakemake.output.plot2, width=800, height=300)
