import pandas as pd
import plotly.graph_objs as go


csv_file_paths = str(snakemake.input).split(" ")
plot_names = ["Barcode", "Cutadapt", "NanoFilt", "Guppy"]

# set a color list for easier viewing
color_list = [
    "rgb(111, 129, 246)",  # barcode
    "rgb(178, 178, 102)",  # cutadapt
    "rgb(61, 194, 156)",   # nanofilt
    "rgb(163, 113, 244)"]  # guppy

# create a list of data frames that match the follow the following list: barcode, cutadapt, nanofilt, guppy
data_frames = []
for file in csv_file_paths:
    data_frames.append(pd.read_csv(
        filepath_or_buffer=file,
        header=0,
        delimiter="\t"
    ))

"""
We are going to iterate through each data frame generated in the list above
For each data frame [frame], we are going to locate [.loc] every instance where the word `unclassified` appears in the column `barcode`
    Then, we will take the corresponding value that is listed under the `reads` column, and extract the value only using `.tolist()[0]`
"""
unclassified_reads_value = []
for frame in data_frames:
    try:
        unclassified_reads_value.append(frame.loc[(frame["barcode"] == "unclassified", "reads")].tolist()[0])
    except IndexError:
        unclassified_reads_value.append("ERROR: No data in .csv file")

# We want to remove the `unclassified` rows from the data frames as they do not let us see the results clearly. Simply remove the last row
for frame in data_frames:
    frame = frame.drop(frame.tail(1).index, inplace=True)

# create a figure to hold our box plots
box_plot = go.Figure()

"""
A varying number of reads can exist outside the 'whiskers' of the plot
We are going to append an appropriate label to the hover_templates_list depending on the number of barcodes present
    i.e. No barcodes, 1 barcode, X barcodes
The hover_templates_list will look as follows:
    [
        [
            "No barcodes",
            "5 barcodes",
            etc.
        ],
        [
            "10 barcodes',
            "1 barcode",
            "No barcodes",
            etc.
        ]
    ]
The end results is a list containing four lists
Each of the inner lists is responsible for one graph (barcode, cutadapt, filtering, mapping)
"""
hover_templates_list = []
for frame in data_frames:
    temp_template = []
    for read in frame["reads"]:
        if read == 0:
            temp_template.append("No barcodes")
        elif read == 1:
            temp_template.append("1 barcode")
        else:
            temp_template.append(f"{read} barcodes")
    hover_templates_list.append(temp_template)


"""
We are now adding 'traces' to the figure
Each trace is an additional box plot; in total there will be four traces (one for each name listed in `plot_names` above
"""
for index, name in enumerate(plot_names):
    box_plot.add_trace(
        go.Box(
            name=name,  #................................................ Set x-axis name
            y=data_frames[index]["reads"],  #............................ Set data for y-axis
            jitter=1.0,  #............................................... This will prevent data points from overlapping (set between 0, 1)
            boxpoints="all",  #.......................................... Show all data points
            marker=dict(  # .............................................. Modify outliers
                color=color_list[index],
                outliercolor="rgb(153, 153, 153)",
                line=dict(
                    outliercolor="rgb(255, 51, 51)",
                    outlierwidth=1)),
            pointpos=0,  #............................................... Set position of data points relative to trace
            hovertext=[label for label in hover_templates_list[index]]  # Set the text for hovering
        ))

# modify appearance of hovering data
hover_label = dict(
    bgcolor="white",
    font_size=16,
    font_family="Rockwell")

# set a title, its location, and size
title_data = dict(
    text="Reads per Barcode after Pipeline Steps",
    x=0.48,
    xanchor="center",
    font=dict(size=25)
)

# create a list of our annotation data; a comment is above each `dict` specifying what the annotation is for
annotation_data = [
    # create a sub title
    dict(xref="paper",
         yref="paper",
         x=0.5,
         y=1.031,
         showarrow=False,
         text="Each count performed after process listed on x-axis",
         font=dict(size=14)),

    # barcode unclassified reads
    dict(xref="paper",
         yref="paper",
         x=0.07,
         y=0.01,
         showarrow=False,
         text=f"Unclassified Reads: {unclassified_reads_value[0]}",
         font=dict(size=13)),

    # cutadapt unclassified reads
    dict(xref="paper",
         yref="paper",
         x=0.375,
         y=0.01,
         showarrow=False,
         text=f"Unclassified Reads: {unclassified_reads_value[1]}",
         font=dict(size=13)),

    # filtering unclassified reads
    dict(xref="paper",
         yref="paper",
         x=0.625,
         y=0.01,
         showarrow=False,
         text=f"Unclassified Reads: {unclassified_reads_value[2]}",
         font=dict(size=13)),

    # mapping unclassified reads
    dict(xref="paper",
         yref="paper",
         x=0.93,
         y=0.01,
         showarrow=False,
         text=f"Unclassified Reads: {unclassified_reads_value[3]}",
         font=dict(size=13)),


]

# update the layout for the graph
box_plot.update_layout(
    hoverlabel=hover_label,
    title=title_data,
    annotations=annotation_data,
    xaxis_title="Process in Pipeline",
    yaxis_title="Reads per Barcode",
    font=dict(size=16)
)

box_plot.write_html(str(snakemake.output))
