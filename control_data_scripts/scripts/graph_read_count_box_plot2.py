'''
Special purpose script for graphing control sequences. Goal is to see how
consistent read counts are between samples within a concentration, and between
concentrations. Concentrations are hard coded as 1, 10, 100, 1000, and 10000.
Read counts are drawn from samp_read_counts.tsv (output by
seekdeep_analysis_snakemake).

Data gets reformatted as a list of concentrations (one item per sample) and a
list of haplotype counts for each amplicon (again, one count per sample).
Iterates through each amplicon and graphs the haplotype count for each sample
within a given concentration.

Graph was drawn with code from here:
https://community.plotly.com/t/how-to-show-overlap-points-in-scatter-plot/24148/13

And here:
https://stackoverflow.com/questions/46750462/subplot-with-plotly-with-multiple-traces

And here:
https://community.plotly.com/t/solved-how-to-create-subplots-using-plotly-express/52418

And here:
https://stackoverflow.com/questions/63460213/how-to-define-colors-in-a-figure-using-plotly-graph-objects-and-plotly-express

import numpy as np
import plotly.graph_objs as go
import plotly

a = np.random.normal(0,1,100)
b = np.random.normal(-2,5,100)

c = np.random.normal(0,1,100)
d = np.random.normal(-2,5,100)

fig = plotly.tools.make_subplots(rows=2,cols=1)

trace_rm1 = go.Histogram(x = a, opacity = 0.75, name = 'malignant')
trace_rm2 = go.Histogram(x = b, opacity = 0.75, name = 'benign')
fig.append_trace(go.Histogram(x = a, opacity = 0.75, name = 'benign'),1,1)
fig.append_trace(go.Histogram(x = b, opacity = 0.75, name = 'malignant'),1,1)
fig.append_trace(go.Histogram(x = c, opacity = 0.75, name = 'benign'),2,1)
fig.append_trace(go.Histogram(x = d, opacity = 0.75, name = 'malignant'),2,1)
fig.layout.update(go.Layout(barmode = 'overlay',))

plotly.offline.plot(fig)

This version uses subplots instead of facets
'''


import plotly.express as px
import pandas as pd
import subprocess
import plotly.subplots as sp
import plotly.graph_objects as go
read_path='/nfs/jbailey5/baileyweb/asimkin/other_people/jmsadler/bigv10ctrl/analyzed_by_PCR_replicates/seekdeep_analysis_snakemake/half_percent_exact_only_analysis_graphs/samp_read_counts.tsv'
output_folder='read_count_box_plots'

desired_concentrations=['1', '10', '100', '1000', '10000']
concentrations, samples=[],[]
conc_dict={}
for line_number, line in enumerate(open(read_path)):
	line=line.strip().split('\t')
	if line_number==0:
		h_dict={}
		for column_number, column in enumerate(line):
			h_dict[column]=column_number
	else:
		conc=line[0].split('-')[0]
		if conc=='10K':
			conc='10000'
		if conc in desired_concentrations:
			samples.append(line[0])
			concentrations.append(conc)
			for column in h_dict:
				if column!='sample':
					conc_dict.setdefault(column, []).append(float(line[h_dict[column]]))

subprocess.call(['mkdir', output_folder])

amplicons=sorted(list(conc_dict.keys()))
fig = sp.make_subplots(rows=6,cols=4, subplot_titles=amplicons, horizontal_spacing = 0.04, vertical_spacing = 0.04)
for amplicon_number, amplicon in enumerate(amplicons):
	data={"concentration": concentrations, "log2_read_count": conc_dict[amplicon], "samples": samples}
	df = pd.DataFrame(data)
	x_coord=amplicon_number//6+1
	y_coord=amplicon_number%6+1
	express_subplot=px.strip(df, x='concentration', y='log2_read_count', hover_data=['samples']).update_traces(jitter = 1).update_traces(marker=dict(color='black'))
	dm = df.groupby('concentration').mean(numeric_only=True)
	ds = df.groupby('concentration').std(numeric_only=True)
	error_subplot=px.scatter(x=dm.index, y=dm['log2_read_count'], 
	error_y=ds['log2_read_count']).update_traces(marker=dict(color='red'))
#	express_subplot.add_scatter(x=dm.index, y=dm['log2_read_count'], 
#	error_y_array=ds['log2_read_count'],
#	mode='markers', showlegend=False).update_traces(marker=dict(color='black'))
	trace_subplot=[]
	print(x_coord, y_coord)
	for trace in range(len(express_subplot["data"])):
		trace_subplot.append(express_subplot["data"][trace])
	for trace in range(len(error_subplot["data"])):
		trace_subplot.append(error_subplot["data"][trace])
	for trace in trace_subplot:
		fig.append_trace(trace, row=y_coord, col=x_coord)
		fig.update_yaxes(title_text="log2 read count", range=[-1, 20], row=y_coord, col=x_coord)
		fig.update_xaxes(title_text="parasitemia (parasites/ul)", row=y_coord, col=x_coord)
#	fig.show()
fig.update_layout(height=3000, width=2000)
fig.write_html(f'{output_folder}/all_strip_plots_subplots.html')
fig.write_image(f'{output_folder}/all_strip_plots_subplots.svg')