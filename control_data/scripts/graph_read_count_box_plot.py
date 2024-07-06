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
for amplicon in conc_dict:
	data={"concentration": concentrations, "log2_read_count": conc_dict[amplicon], "samples": samples}
	df = pd.DataFrame(data)
	fig=px.box(df, x="concentration", y="log2_read_count", points="all", hover_data=['samples'])
	fig['layout']['title']=amplicon
#	print(fig)
#	print(amplicon)
#	fig.show()
	fig.write_html(f'{output_folder}/{amplicon}_box_plot.html')
#	exit()