'''
reports the number of amplicons successfully sequenced for every sample
'''
import yaml
import plotly.express as px
import pandas as pd

yaml_file='/nfs/jbailey5/baileyweb/asimkin/other_people/jmsadler/bigv10ctrl/analyzed_by_PCR_replicates/seekdeep_analysis_snakemake/half_percent_exact_only_analysis_graphs/read_dict.yaml'
main_db=yaml.safe_load(open(yaml_file))
output_folder='amplicons_per_sample'

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
	fig=px.box(df, x="parasitemias", y="amplicons retrieved", points="all", hover_data=['samples'])
	fig['layout']['title']='amplicons retrieved per sample at each parasitemia'
	fig.write_html(f'{output_folder}/amplicons_per_sample_box_plot.html')
#	exit()

reorganized_db=reorganize_db(main_db)
graph_db(reorganized_db)