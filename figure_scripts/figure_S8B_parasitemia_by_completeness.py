'''
This uses Jacob's parasitemia spreadsheet plus my tabular output of S8A to
create S8B. Jacob's parasitemia spreadsheet is from here:
https://docs.google.com/spreadsheets/d/1BSyIwm7JbqAniuAfz36fX7m14QlFl6lu/edit?
gid=873760498#gid=873760498
'''
output_folder='better_log_scale'
parasitemia_values=('/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_'
	'admin/pfsmarrt_github/field_data/Full_data_Cameroon-DSG2020_Speciation_'
	'Summary.tsv')

completeness_values=('/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_'
	'admin/pfsmarrt_github/field_data/manuscript_figures/S8A_amplicons_per_'
	'sample.tsv')

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

def graph_db():
	'''
	flattens the database and reformats it into a somewhat redundant format for
	graphing in plotly
	'''
	import pandas as pd
	import math
	import plotly.express as px
	import subprocess
	import statistics
	subprocess.call(['mkdir', output_folder])
	sorted_samples=sorted(list(parasitemia_dict.keys()))
	y_values, log_y_values, categories=[],[],[]
	for sample in sorted_samples:
		log_y_values.append(math.log(float(parasitemia_dict[sample]), 10))
		y_values.append(float(parasitemia_dict[sample]))
		categories.append(completeness_dict[sample])
	data={'categories': categories, 'samples': sorted_samples, 'log10 parasitemia': log_y_values}
	sd_dict={}
	for category_number, category in enumerate(categories):
		sd_dict.setdefault(category, []).append(y_values[category_number])
	cat_list, mean_list, stdev_list=[],[],[]
	for category in sd_dict:
		mean=round(statistics.mean(sd_dict[category]), 3)
		stdev=round(statistics.stdev(sd_dict[category]), 3)
		print(category, mean, stdev)
		mean_list.append(math.log(mean, 10))
		cat_list.append(category)
		stdev_list.append(math.log(stdev, 10))
	write_table(data, f'{output_folder}/S8B_parasitemia_by_completeness.tsv')
	df = pd.DataFrame(data)
#	fig=px.strip(df, hover_data=['samples']).update_traces(jitter = 1).update_traces(marker=dict(color='black'))
	fig=px.strip(df, x='categories', y='log10 parasitemia', hover_data=['samples']).update_traces(marker=dict(color='black'))
	#dm = df.groupby('categories').mean(numeric_only=True)
	#ds = df.groupby('categories').std(numeric_only=True)
	#need to manually take log of non-log means and standard deviations here
	#print('dm is', dm)
	#print('ds is', ds)

	#use plotly express to create a scatterplot of error bars, and update all
	#traces to be red
	
	fig.add_scatter(x=cat_list, y=mean_list, 
		error_y_array=stdev_list,
        	mode='markers', showlegend=False, fillcolor='red')
	#fig.update_traces(marker=dict(color='red'))
	fig.write_html(f'{output_folder}/S8B_parasitemia_by_completeness.html')
	fig.write_image(f'{output_folder}/S8B_parasitemia_by_completeness.svg')
	fig.write_image(f'{output_folder}/S8B_parasitemia_by_completeness.png')
	return df


parasitemia_dict=dict([line.strip().split('\t') for line in open(parasitemia_values)][1:])
completeness_list=[line.strip().split('\t') for line in open(completeness_values)][1:]
completeness_dict=dict([[line[-1].replace('-', '_'), line[0]] for line in completeness_list])
for sample in completeness_dict:
	if completeness_dict[sample]=='24':
		completeness_dict[sample]='complete'
	else:
		completeness_dict[sample]='incomplete'
graph_db()