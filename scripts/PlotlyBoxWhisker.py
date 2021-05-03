import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots

csv_file_paths = str(snakemake.input).split(" ")
plot_names = []  # only generate plots we have data for
color_list = []  # set a color list for easier viewing

for i, file in enumerate(csv_file_paths):
    if "barcode" in file:
        plot_names.append("Barcode (Guppy Barcoder)")
        color_list.append("rgb(111, 129, 246)")
    elif "cutadapt" in file:
        plot_names.append("Trim (Cutadapt)")
        color_list.append("rgb(178, 178, 102)")
    elif "filter" in file:
        plot_names.append("Filtering (NanoFilt)")
        color_list.append("rgb(61, 194, 156)")
    elif "mapping" in file:
        plot_names.append("Mapping (MiniMap)")
        color_list.append("rgb(163, 113, 244)")

# create a list of data frames for our input files
data_frames = []
for file in csv_file_paths:
    data_frames.append(pd.read_csv(
        filepath_or_buffer=file,
        header=0,
        delimiter="\t"
    ))

"""
We are going to iterate through each data frame in the data_frames list above
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

# Create a subplot figure to hold our graphs
box_plot = make_subplots(rows=1,
                         cols=len(plot_names),
                         x_title="Process in Pipeline",
                         y_title="Reads per Barcode",
                         shared_yaxes=True,
                         shared_xaxes=True,
                         horizontal_spacing=0)

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
Each trace is an additional box plot; in total there will be four traces (or one for each name listed in `plot_names` above, if not all graphs available)
"""
for index, name in enumerate(plot_names):
    box_plot.add_trace(
        go.Box(
            name=name,  # ................................................ Set x-axis name
            y=data_frames[index]["reads"],  # ............................ Set data for y-axis
            jitter=1.0,  # ............................................... This will prevent data points from overlapping (set between 0, 1)
            boxpoints="all",  # .......................................... Show all data points
            marker=dict(  # .............................................. Modify outliers
                color=color_list[index],
                outliercolor="rgb(153, 153, 153)",
                line=dict(
                    outliercolor="rgb(255, 51, 51)",
                    outlierwidth=1)),
            pointpos=0,  # ............................................... Set position of data points relative to trace
            hovertext=[label for label in hover_templates_list[index]]  # Set the text for hovering
        ),
        row=1, col=index+1
    )

# remove legend from plot
for i in range(len(plot_names)):
    box_plot.update_traces(row=1, col=i+1, showlegend=False)

# add annotations
# Subtitle annotation
box_plot.add_annotation(
    dict(
        xref="paper",
        yref="paper",
        x=0.,
        y=1.045,
        showarrow=False,
        text="Each count performed after process listed on x-axis",
        font=dict(size=14)
    )
)

# Unclassified reads annotations
for i, read in enumerate(unclassified_reads_value):
    box_plot.add_annotation(
        dict(
            xref=f"x{i+1}",
            yref=f"y{i+1}",
            x=0,
            y=1,
            text=f"Unclassified Reads: {read}",
            showarrow=False
        )
    )

# update the layout for the graph
box_plot.update_layout(
    # set mouse-hovering data
    hoverlabel=dict(
        bgcolor="white",
        font_size=16,
        font_family="Rockwell"),
    # set title datas
    title=dict(
        text="Reads per Barcode after Pipeline Steps",
        xanchor="auto",
        font=dict(size=25)),
    font=dict(size=16)
)

box_plot.update_yaxes(automargin=True)
box_plot.update_xaxes(automargin=True)

box_plot.write_html(str(snakemake.output))
