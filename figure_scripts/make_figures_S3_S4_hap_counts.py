'''
Special purpose script for graphing control sequences. Goal is to see if
different concentrations of 3D7 always yield one haplotype in seekdeep.
Concentrations are hard coded as 1, 10, 100, 1000, and 10000. Haplotype counts
are drawn from samp_hap_counts.tsv (output by seekdeep_analysis_snakemake).

Data gets reformatted as a list of concentrations (one item per sample) and a
list of haplotype counts for each amplicon (again, one count per sample).
Iterates through each amplicon and graphs the haplotype count for each sample
within a given concentration.
'''

import plotly.express as px
import pandas as pd
import subprocess
import plotly

COI_path='/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_admin/seekdeep_graphs/controls_0.5_output/samp_hap_counts.tsv'
output_folder='manuscript_figures'

desired_concentrations=['1', '10', '100', '1000', '10000']
concentrations, samples=[],[]
conc_dict={}
for line_number, line in enumerate(open(COI_path)):
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

amplicons=sorted(list(conc_dict.keys()))
subprocess.call(['mkdir', output_folder])


def graph_amplicons(amplicon_list, grid_rows, grid_columns, x_dimensions, y_dimensions, plot_name):
	#create empty subplot figure
	fig = plotly.subplots.make_subplots(rows=grid_rows,cols=grid_columns, subplot_titles=amplicon_list, horizontal_spacing = 0.04, vertical_spacing = 0.06)
	for amplicon_number, amplicon in enumerate(amplicon_list):
		data={"concentration": concentrations, "haplotype count": conc_dict[amplicon], "samples": samples}
		df = pd.DataFrame(data)
		x_coord=amplicon_number%grid_columns+1
		y_coord=amplicon_number//grid_columns+1

		#use plotly express to make a strip plot of values at each concentration,
		#and update all traces in it to be black. Also maximize jitter for easier
		#point differentiation
		express_subplot=px.strip(df, x='concentration', y='haplotype count', hover_data=['samples']).update_traces(jitter = 1).update_traces(marker=dict(color='black'))

		#add all the traces from the strip plot to a list
		trace_subplot=[]
		print(x_coord, y_coord)
		for trace in range(len(express_subplot["data"])):
			trace_subplot.append(express_subplot["data"][trace])

		#append all traces onto the relevant square of the subplot figure
		for trace in trace_subplot:
			fig.append_trace(trace, row=y_coord, col=x_coord)

		#update the axes of the relevant subplot figure
		fig.update_yaxes(title_text="haplotype count", range=[-0.3, 4], row=y_coord, col=x_coord)
		fig.update_xaxes(title_text="parasitemia (parasites/ul)", row=y_coord, col=x_coord)

	#update the layout so the graph is appropriately sized without weird stretching
	fig.update_layout(height=y_dimensions, width=x_dimensions)
	fig.write_html(f'{output_folder}/{plot_name}.html')

	#write_image requires kaleido package, but do not use newest version of numpy
	# (breaks the program for some reason - I reverted to 1.26.4) 
	fig.write_image(f'{output_folder}/{plot_name}.svg')
	fig.write_image(f'{output_folder}/{plot_name}.png')

heome_amplicons=['ama1', 'heome-a', 'heome-b', 'heome-c', 'heome-d', 'heome-e',
'heome-f', 'heome-g', 'heome-h']
graph_amplicons(heome_amplicons, 3, 3, 1500, 1500, 'S3_heome_haps')

print('moving on')

DR_amplicons=['dhfr-108', 'dhfr-51-59', 'dhps-436-437', 'dhps-540', 'dhps-581',
'dhps-613', 'k13-a', 'k13-b', 'k13-c', 'k13-f', 'k13-g', 'mdr1-1034',
'mdr1-184', 'mdr1-86', 'pfcrt']
graph_amplicons(DR_amplicons, 5, 3, 1500, 2500, 'S4_DR_haps')