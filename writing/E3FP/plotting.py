#!/usr/bin/env python3
import numpy as np
import pandas as pd
import plotly.express as px

distances = np.loadtxt("E3FP.mat")
projection = pd.read_table("MDS3.tsv")

# yellow-green-blue color scale that fits with the other colors in the fig
heatmap = px.imshow(distances, color_continuous_scale=px.colors.sequential.Viridis)
heatmap.update_xaxes(showticklabels=False)
heatmap.update_yaxes(showticklabels=False)
heatmap.show()
heatmap.write_image("E3FP_heatmap.svg", width=3, height=5)

labels = dict(MDS_1="MDS 1", MDS_2="MDS 2", MDS_3="MDS 3")
proj = px.scatter_3d(projection, x="MDS_1", y="MDS_2", z="MDS_3", labels=labels, hover_name="id", width=600, height=600)
proj.update_traces(marker_size=4)
# write_image doesn't zoom well with 3D


