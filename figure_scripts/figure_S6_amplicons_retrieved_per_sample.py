'''
reports the number of amplicons successfully sequenced for every sample
'''
import yaml
import plotly.express as px
import pandas as pd

yaml_file='/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_admin/seekdeep_graphs/controls_0.5_output/read_dict.yaml'
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
			conc=sample.split('-')[0]
			if conc in ['1', '10', '100', '1000', '10K']:
				if conc=='1000':
					conc='1,000'
				if conc=='10K':
					conc='10,000'
				reorganized_db.setdefault(conc, {})
				new_sample=sample.replace('-', '.')
				found_amplicon=False
				for replicate in main_db[amplicon][sample]:
					if len(main_db[amplicon][sample][replicate])>0:
						found_amplicon=True
				if found_amplicon:
					reorganized_db[conc][new_sample]=reorganized_db[conc].setdefault(new_sample, 0)+1
	return reorganized_db

def graph_db(reorganized_db):
	'''
	flattens the database and reformats it into a somewhat redundant format for
	graphing in plotly
	'''
	import subprocess
	subprocess.call(['mkdir', output_folder])
	concentrations, samples, counts=[],[],[]
	for conc in reorganized_db:
		for sample in reorganized_db[conc]:
			concentrations.append(conc)
			samples.append(sample)
			counts.append(reorganized_db[conc][sample])
	data={'samples': samples, 'amplicons retrieved': counts, 'parasitemias': concentrations}
	df = pd.DataFrame(data)

	fig=px.strip(df, x='parasitemias', y='amplicons retrieved', hover_data=['samples']).update_traces(jitter = 1).update_traces(marker=dict(color='black'))

	dm = df.groupby('parasitemias').mean(numeric_only=True)
	ds = df.groupby('parasitemias').std(numeric_only=True)

	#use plotly express to create a scatterplot of error bars, and update all
	#traces to be red
	fig.add_scatter(x=dm.index, y=dm['amplicons retrieved'], 
		error_y_array=ds['amplicons retrieved'],
        	mode='markers', showlegend=False, fillcolor='red')
	#fig.update_traces(marker=dict(color='red'))
	fig.update_xaxes(title_text="parasitemia (parasites/ul)")
	fig.write_html(f'{output_folder}/S6_amplicons_per_sample.html')
	fig.write_image(f'{output_folder}/S6_amplicons_per_sample.svg')
	fig.write_image(f'{output_folder}/S6_amplicons_per_sample.png')

#	exit()

reorganized_db=reorganize_db(main_db)
graph_db(reorganized_db)