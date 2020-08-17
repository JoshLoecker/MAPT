import numpy as np
import pandas as pd
import plotly.express as px

csv_file = snakemake.input[0]

save_file = snakemake.output[0]

# load our csv file
data_frame = pd.read_csv(
    filepath_or_buffer=csv_file,
    delimiter="\t",
    header=0)

# we want to collect the number of unclassified reads to create an annotation on the graph
try:
    unclassified_reads = data_frame.loc[ data_frame.barcode == "unclassified",'reads' ].tolist()[0]
except IndexError:
    unclassified_reads = "ERROR: No data in .csv file"

# remove the unclassified row from the frame. It is much larger than anything else and does not allow us to see data properly
data_frame_remove_unclassified = data_frame
data_frame_remove_unclassified.drop( data_frame_remove_unclassified.tail(1).index, inplace=True )

# get the average reads with and without the unclassified reads
try:
    average_reads_with_unclassified = int(np.average(data_frame['reads']))
    average_reads_without_classified = int(np.average(data_frame_remove_unclassified['reads']))
except (ZeroDivisionError, ValueError):
    average_reads_with_unclassified = "ERROR: No data in .csv file"
    average_reads_without_unclassified = "ERROR: No data in .csv file"


# create the histogram
histogram = px.histogram(data_frame=data_frame_remove_unclassified,
                         x="barcode",
                         y="reads",
                         marginal="box",
                         hover_name=data_frame_remove_unclassified['barcode'])

# modify the label when hovering over a bar
histogram.update_traces(hovertemplate="%{x}<br>%{y} reads")

# modify data on hovering over chart
hover_label = dict(
    bgcolor="white",
    font_size=16,
    font_family="Rockwell"
)

# set a title, its location, and size
title_data = dict(
    text="Reads per Barcode",
    x=0.5,
    xanchor="center",
    font=dict(size=25)
)

# create a list of our annotation data
# a comment is above each `dict` specifying what the annotation is for
annotation_data = [
    # create a sub title
    dict(xref="paper",
         yref="paper",
         x=0.5,
         y=1.031,
         showarrow=False,
         text=str(snakemake.params.sub_title),
         font=dict(size=18)),

    # show the number of unclassified reads
    dict(xref="paper",
         yref="paper",
         x=1.00,
         y=1.03,
         showarrow=False,
         text="Unclassified Reads: {0}".format(unclassified_reads),
         font=dict(size=18)),

    # show the average number of reads, excluding unclassified
    dict(xref="paper",
         yref="paper",
         x=1.00,
         y=1.00,
         showarrow=False,
         text="Average reads (excluding unclassified): {0}".format(average_reads_without_unclassified),
         font=dict(size=18))
]


# update the layout for the graph
histogram.update_layout(
    hoverlabel=hover_label,
    title=title_data,
    annotations=annotation_data,
    xaxis_title="Barcode Number",
    yaxis_title="Reads per Barcode",
    font=dict(size=14)
)

histogram.write_html(save_file)
