'''
stacked bar charts are not supported in plotly graph objects (as far as I can
tell) and plotly express (stacked bar charts) don't seem to be supported by
subplots (which want graph object style traces). It might be possible (with
plotly graph objects) to fool the system into plotting bars on top of each
other. This website looks promising:
https://stackoverflow.com/questions/70563166/stacked-barplot-in-plotly
(has graph object bar charts that are stacked)
'''
import yaml
import plotly.express as px
import pandas as pd
import plotly
import plotly.graph_objects as go

yaml_file='/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_admin/seekdeep_graphs/controls_0.5_output/read_dict.yaml'
main_db=yaml.safe_load(open(yaml_file))

def reorganize_db(main_db):
	'''
	sums haplotype counts across replicates within each sample, keeps only
	samples that are numeric concentrations
	'''
	reorganized_db={}
	for amplicon in main_db:
		for sample in main_db[amplicon]:
			if sample.split('-')[0] in ['1', '10', '100', '1000', '10K']:
				new_sample='s_'+sample.replace('-', '.')
				hap_dict={}
				for replicate in main_db[amplicon][sample]:
					for hap in main_db[amplicon][sample][replicate]:
						new_hap=hap.split('.')[1].replace('4', '1')
						hap_count=main_db[amplicon][sample][replicate][hap]
						hap_dict[new_hap]=hap_dict.setdefault(new_hap, 0)+hap_count
				if amplicon not in reorganized_db:
					reorganized_db[amplicon]={}
				reorganized_db[amplicon][new_sample]=hap_dict
	return reorganized_db

def get_fracs(hap_counts):
	'''
	converts counts for each haplotype into fractions
	'''
	frac_dict={}
	for amplicon in hap_counts:
		for sample in hap_counts[amplicon]:
			hap_total=0
			for hap in hap_counts[amplicon][sample]:
				hap_total+=hap_counts[amplicon][sample][hap]
			for hap in hap_counts[amplicon][sample]:
				hap_frac=hap_counts[amplicon][sample][hap]/hap_total
				#graph only 'unexpected' haplotypes
				if hap_frac<1 or hap.split('.')[-1]!='0':
					frac_dict.setdefault(amplicon, {})
					frac_dict[amplicon].setdefault(sample, {})
					frac_dict[amplicon][sample][hap]=hap_frac
	return frac_dict

def flatten_graphing(frac_counts):
	'''
	creates a repetitive data structure suitable for graphing, where every
	fraction has a redundantly labeled sample and haplotype
	'''
	graphing_dict={}
	for amplicon in frac_counts:
		samples, haps, fracs=[],[],[]
		for sample in frac_counts[amplicon]:
			for hap in frac_counts[amplicon][sample]:
				samples.append(sample)
				haps.append(hap)
				fracs.append(frac_counts[amplicon][sample][hap])
		data={'samples': samples, 'haps': haps, 'fracs': fracs}
		df = pd.DataFrame(data)
		graphing_dict[amplicon]=df
	return graphing_dict

def graph_plots(amplicon_list, grid_rows, grid_columns, x_dimensions, y_dimensions, plot_name):
	import subprocess
	subprocess.call(['mkdir', 'stacked_bars'])
	fig = plotly.subplots.make_subplots(rows=grid_rows,cols=grid_columns, subplot_titles=amplicons, horizontal_spacing = 0.04, vertical_spacing = 0.30)
	for amplicon_number, amplicon in enumerate(amplicon_list):
		x_coord=amplicon_number%grid_columns+1
		y_coord=amplicon_number//grid_columns+1
		df=graphing_dict[amplicon]
		express_subplot=px.bar(df, x="samples", y="fracs", color="haps", title=amplicon).update_traces(width = 0.1)
#		print('express subplot is', express_subplot)
		#add all the traces from the strip plot to a list
		trace_subplot=[]
		print(x_coord, y_coord)
		for trace in range(len(express_subplot["data"])):
			trace_subplot.append(express_subplot["data"][trace])


		#append all traces onto the relevant square of the subplot figure
		for trace in trace_subplot:
			fig.append_trace(trace, row=y_coord, col=x_coord)

		#update the axes of the relevant subplot figure
		fig.update_yaxes(title_text="haplotype fraction", row=y_coord, col=x_coord)
		fig.update_xaxes(title_text="sample name", row=y_coord, col=x_coord)
		fig.update_traces(offsetgroup='0')
		fig.update_traces(legendgroup='0')
#		fig.update_traces(base=express_subplot['data'][0])
		print('fig is', fig)
	fig.write_html('stacked_bars/'+plot_name+'.html')
	fig.write_image('stacked_bars/'+plot_name+'.svg')
	fig.write_image('stacked_bars/'+plot_name+'.png')

hap_counts=reorganize_db(main_db)
frac_counts=get_fracs(hap_counts)
graphing_dict=flatten_graphing(frac_counts)
amplicons=sorted(list(graphing_dict.keys()))

#print('ama1 is', graphing_dict['ama1'])
#print('medals is', px.data.medals_long())

#medals = px.bar(px.data.medals_long(), x="nation", y="count", color="medal", title="Long-Form Input")
#print('properly stacked is', medals)
#medals.write_html('medals.html')

print('keys are', graphing_dict.keys())

for amplicon in graphing_dict:
	haps=px.bar(graphing_dict[amplicon], x="samples", y="fracs", color="haps", title=amplicon).update_traces(width = 0.3)
	print('haps is', haps)
	haps.write_html(f'{amplicon}_fracs.html')

graph_plots(amplicons, 2, 2, 1000, 1000, 'all_unexpected_hap_fracs')

