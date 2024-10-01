'''
reports the number of amplicons successfully sequenced for every sample
'''
import yaml
import plotly.express as px
import pandas as pd

yaml_file='/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_admin/seekdeep_graphs/ind_field_samples_new_replicates_new_params_389/read_dict.yaml'
main_db=yaml.safe_load(open(yaml_file))
output_folder='manuscript_figures'

def reorganize_db(main_db):
	'''
	sums haplotype counts across replicates within each sample, keeps only
	samples that are numeric concentrations
	'''
	reorganized_db={}
	for amplicon in main_db:
		for sample in main_db[amplicon]:
			found_amplicon=False
			for replicate in main_db[amplicon][sample]:
				if len(main_db[amplicon][sample][replicate])>0:
					found_amplicon=True
			if found_amplicon:
				reorganized_db[sample]=reorganized_db.setdefault(sample, 0)+1
	return reorganized_db

def write_table(data, output_path):
	output_file=open(output_path, 'w')
	output_line=[]
	keys=sorted(list(data.keys()))
	for key in keys:
		output_line.append(key)
	output_file.write('\t'.join(output_line)+'\n')
	for count_number in range(len(data[keys[0]])):
		output_line=[]
		for key in keys:
			output_line.append(str((data[key][count_number])))
		output_file.write('\t'.join(output_line)+'\n')

def graph_db(reorganized_db):
	'''
	flattens the database and reformats it into a somewhat redundant format for
	graphing in plotly
	'''
	import subprocess
	subprocess.call(['mkdir', output_folder])
	titles, counts, samples=[],[],[]
	for sample in reorganized_db:
		counts.append(reorganized_db[sample])
		samples.append(sample)
		titles.append('all field samples')
	data={'amplicons retrieved per sample': counts, 'samples':samples, 'category':titles}
	write_table(data, f'{output_folder}/S8A_amplicons_per_sample.tsv')
	df = pd.DataFrame(data)
#	fig=px.strip(df, hover_data=['samples']).update_traces(jitter = 1).update_traces(marker=dict(color='black'))
	fig=px.strip(df, x='category', y='amplicons retrieved per sample', hover_data=['samples']).update_traces(marker=dict(color='black'))

	dm = df.groupby('category').mean(numeric_only=True)
	ds = df.groupby('category').std(numeric_only=True)

	#use plotly express to create a scatterplot of error bars, and update all
	#traces to be red
	fig.add_scatter(x=dm.index, y=dm['amplicons retrieved per sample'], 
		error_y_array=ds['amplicons retrieved per sample'],
        	mode='markers', showlegend=False, fillcolor='red')
	#fig.update_traces(marker=dict(color='red'))
	fig.update_xaxes(title_text="amplicons retrieved per sample")
	fig.write_html(f'{output_folder}/S8A_amplicons_per_sample.html')
	fig.write_image(f'{output_folder}/S8A_amplicons_per_sample.svg')
	fig.write_image(f'{output_folder}/S8A_amplicons_per_sample.png')
	return df

reorganized_db=reorganize_db(main_db)
data=graph_db(reorganized_db)