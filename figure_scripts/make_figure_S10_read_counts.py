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
read_path='/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_admin/seekdeep_graphs/controls_0.5_output/samp_read_counts.tsv'
output_folder='manuscript_figures'

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


#create empty subplot figure
fig = sp.make_subplots(rows=6,cols=4, subplot_titles=amplicons, horizontal_spacing = 0.04, vertical_spacing = 0.04)

#each amplicon has the same concentrations and the same samples (5
#concentrations with 12 samples each, for 60 concentration, sample pairs
#(represented as two lists of 60 elements each). Approach here is to combine
#these two lists with the 60 unique points that represent each amplicon, to form
#a pandas dataframe of 3 lists with 60 points each
for amplicon_number, amplicon in enumerate(amplicons):
	data={"concentration": concentrations, "log2_read_count": conc_dict[amplicon], "samples": samples}
	df = pd.DataFrame(data)
	x_coord=amplicon_number%4+1
	y_coord=amplicon_number//4+1
	print('x is', x_coord, 'y is', y_coord)
	#use plotly express to make a strip plot of values at each concentration,
	#and update all traces in it to be black. Also maximize jitter for easier
	#point differentiation
	express_subplot=px.strip(df, x='concentration', y='log2_read_count', hover_data=['samples']).update_traces(jitter = 1).update_traces(marker=dict(color='black'))

	#calculate means and standard deviations for each concentration
	dm = df.groupby('concentration').mean(numeric_only=True)
	ds = df.groupby('concentration').std(numeric_only=True)

	#use plotly express to create a scatterplot of error bars, and update all
	#traces to be red
	error_subplot=px.scatter(x=dm.index, y=dm['log2_read_count'], 
	error_y=ds['log2_read_count']).update_traces(marker=dict(color='red'))

	#add all the traces from the strip plot to a list
	trace_subplot=[]
	for trace in range(len(express_subplot["data"])):
		trace_subplot.append(express_subplot["data"][trace])

	#add all the traces from the error bar plot to the same list
	for trace in range(len(error_subplot["data"])):
		trace_subplot.append(error_subplot["data"][trace])

	#append all traces onto the relevant square of the subplot figure
	for trace in trace_subplot:
		fig.append_trace(trace, row=y_coord, col=x_coord)

	#update the axes of the relevant subplot figure
	fig.update_yaxes(title_text="log2 read count", range=[-1, 20], row=y_coord, col=x_coord)
	fig.update_xaxes(title_text="parasitemia (parasites/ul)", row=y_coord, col=x_coord)

#update the layout so the graph is appropriately sized without weird stretching
fig.update_layout(height=3000, width=2000)
fig.write_html(f'{output_folder}/S2_sequencing_depth.html')

#write_image requires kaleido package, but do not use newest version of numpy
# (breaks the program for some reason - I reverted to 1.26.4) 
fig.write_image(f'{output_folder}/S2_sequencing_depth.svg')
fig.write_image(f'{output_folder}/S2_sequencing_depth.png')