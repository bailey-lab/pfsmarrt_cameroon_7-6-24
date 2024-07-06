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

This version uses facets instead of subplots
'''


import plotly.express as px
import pandas as pd
import subprocess
import plotly
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

counts, amplicons=[],[]
for amplicon in conc_dict:
	counts.extend(conc_dict[amplicon])
	amplicons.extend([amplicon]*len(conc_dict[amplicon]))
concentrations=concentrations*len(conc_dict)
samples=samples*len(conc_dict)

data={"concentration": concentrations, "log2_read_count": counts, "samples": samples, 'amplicon': amplicons}
df = pd.DataFrame(data)
fig=px.strip(df, x='concentration', y='log2_read_count', hover_data=['samples'], facet_col='amplicon', range_x=[0,1000], facet_col_wrap=4, width=4000, height=6000).update_traces(jitter = 1)
#	fig.layout.update(go.Layout(barmode = 'overlay',))
#	dm = df.groupby('concentration').mean()
#	ds = df.groupby('concentration').std()
#	fig.append_trace(px.scatter(x=dm.index, y=dm['log2_read_count'], 
#	error_y_array=ds['log2_read_count'],
#	mode='markers', showlegend=False)
#	fig.show()
fig.write_html(f'{output_folder}/v3_strip_plots.html')