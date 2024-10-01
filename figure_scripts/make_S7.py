'''
converts the parsed outputs of parse_blat.py into a figure S7 representation.
'''
parsed_input='/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_admin/pfsmarrt_github/control_data/map_haplotypes/testing_pipeline/all_haplotypes.tsv'

valid_amps=set(['ama1', 'dhfr-108', 'dhfr-51-59', 'dhps-436-437', 'dhps-613',
	'heome-b', 'heome-c', 'heome-e', 'heome-f', 'heome-g', 'mdr1-86', 'pfcrt'])

def get_haps(amplicon, entry_string):
	entries=entry_string.split(' *** ')
	haps=[entry.split(' ')[0] for entry in entries]
	haps=[hap.replace('no', f'{amplicon}.1') for hap in haps] #this only works because I checked and in this dataset every 'missing' haplotype of a given replicate is a '.1' in other replicates
	fracs=[entry.split(' ')[-1][1:-2] for entry in entries]
	fracs=[float(frac.replace('at', '0')) for frac in fracs] #every 'missing' haplotype is imputed to have a value of 0 for the '.1' version of the haplotype
	return {hap:fracs[hap_number] for hap_number, hap in enumerate(haps)}

#td7 is 3D7, sg8 is 7G8
count_dict={}
data=sorted([line.strip().split('\t') for line in open(parsed_input)][1:])
for line in data:
	replicate, amplicon, hap_count, td7, dd2, sg8=line
	if '955' in replicate and amplicon in valid_amps:
		split_rep=replicate.split('-')
		sample='-'.join(split_rep[:-1])
		td7_dict=get_haps(amplicon, td7)
		dd2_dict=get_haps(amplicon, dd2)
		sg8_dict=get_haps(amplicon, sg8)
		for hap in dd2_dict:
			if hap not in td7_dict:
				count_dict.setdefault(sample, {})
				count_dict[sample].setdefault(hap, []).append(dd2_dict[hap])
#			else:
#				print('hap', hap, 'of rep', replicate, 'is in 3D7 and dd2, with'
#					' frequencies of', td7_dict[hap], 'and', dd2_dict[hap],
#					'respectively')
def graph_data(count_dict):
	import subprocess
	import pandas as pd
	import plotly.express as px
	output_folder='manuscript_figures'
	subprocess.call(f'mkdir {output_folder}', shell=True)
	graphing_table=[]
	for sample in count_dict:
		mdr_summed=0
		for hap in sorted(count_dict[sample]):
			averaged_count=round(sum(count_dict[sample][hap])/len(count_dict[sample][hap]), 2)
			if 'mdr1-86' in hap:
				mdr_summed+=averaged_count
			graphing_table.append([hap, sample, str(averaged_count)])
		graphing_table.append(['mdr1-86_summed', sample, str(mdr_summed)])
	graphing_table.sort()
	x_list=[line[0] for line in graphing_table]
	sample_list=[line[1] for line in graphing_table]
	y_list=[float(line[2]) for line in graphing_table]
	output_file=open(f'{output_folder}/S7_tabular.tsv', 'w')
	for line in graphing_table:
		output_file.write('\t'.join(line)+'\n')
	data={'haplotypes': x_list, 'minor_fractions': y_list, 'samples': sample_list}
	df = pd.DataFrame(data)
	fig=px.strip(df, x='haplotypes', y='minor_fractions', hover_data=['samples']).update_traces(jitter = 1).update_traces(marker=dict(color='black'))
	dm = df.groupby('haplotypes').mean(numeric_only=True)
	ds = df.groupby('haplotypes').std(numeric_only=True)

	#use plotly express to create a scatterplot of error bars, and update all
	#traces to be red
	fig.add_scatter(x=dm.index, y=dm['minor_fractions'], 
		error_y_array=ds['minor_fractions'],
        	mode='markers', showlegend=False, fillcolor='red')
	#fig.update_traces(marker=dict(color='red'))
	fig.update_xaxes(title_text="minor haplotypes")
	fig.update_yaxes(title_text="Dd2 minor fraction (percent)")
	fig.write_html(f'{output_folder}/S7_minor_haplotype_frequencies.html')
	fig.write_image(f'{output_folder}/S7_minor_haplotype_frequencies.svg')
	fig.write_image(f'{output_folder}/S7_minor_haplotype_frequencies.png')
graph_data(count_dict)