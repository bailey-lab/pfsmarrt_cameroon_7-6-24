'''
creates figure S11 from samp_hap_counts.tsv
'''

import plotly.express as px
import pandas as pd
import subprocess
import plotly.subplots as sp
import plotly.graph_objects as go

COI_path='/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_admin/seekdeep_graphs/ind_field_samples_new_replicates_new_params_389/samp_hap_counts.tsv'
output_folder='manuscript_figures'
graphing_path=f'{output_folder}/S11_graphing_table.tsv'

subprocess.call(f'mkdir {output_folder}', shell=True)

def make_graphing_table(COI_path):
	graphing_table=[]
	for line_number, line in enumerate(open(COI_path)):
		line=line.strip().split('\t')
		if line_number==0:
			h_rev_dict={}
			for column_number, column in enumerate(line):
				h_rev_dict[column_number]=column
		else:
			sample=line[0]
			for column_number, hap_count in enumerate(line):
				if column_number>0:
					amplicon=h_rev_dict[column_number]
					graphing_table.append([amplicon, sample, hap_count])
	graphing_table.sort()
	graphing_file=open(graphing_path, 'w')
	graphing_file.write(f'amplicons\tsamples\tcounts\n')
	for line in graphing_table:
		graphing_file.write('\t'.join(list(map(str, line)))+'\n')
	graphing_file.close()

def graph_data(graphing_path):
	import pandas as pd
	import plotly.express as px
	df=pd.read_csv(graphing_path, sep='\t')
	output_folder='manuscript_figures'
	fig=px.strip(df, x='amplicons', y='counts', hover_data=['samples']).update_traces(jitter = 1).update_traces(marker=dict(color='black'))
	dm = df.groupby('amplicons').mean(numeric_only=True)
	ds = df.groupby('amplicons').std(numeric_only=True)
	#use plotly express to create a scatterplot of error bars, and update all
	#traces to be red
	fig.add_scatter(x=dm.index, y=dm['counts'], 
		error_y_array=ds['counts'],
        	mode='markers', showlegend=False, fillcolor='red')
	#fig.update_traces(marker=dict(color='red'))
	fig.update_xaxes(title_text="amplicons")
	fig.update_yaxes(title_text="sample haplotype counts")
	fig.write_html(f'{output_folder}/S11_field_sample_hap_counts.html')
	fig.write_image(f'{output_folder}/S11_field_sample_hap_counts.svg')
	fig.write_image(f'{output_folder}/S11_field_sample_hap_counts.png')
make_graphing_table(COI_path)
graph_data(graphing_path)