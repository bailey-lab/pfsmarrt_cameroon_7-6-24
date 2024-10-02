import subprocess
mccoil_100_smarrt='/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_admin/pfsmarrt_github/abebe_COI_data/Cameroon_PFSMARRTCOI_realmccoi.csv'
mccoil_100_smarrt=[line.strip().split(',') for line in open(mccoil_100_smarrt)][1:]
mccoil_100_smarrt={line[0].replace('_', '-'):float(line[2]) for line in mccoil_100_smarrt}

mccoil_50_mip='/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_admin/pfsmarrt_github/abebe_COI_data/CameroonMIPCOI_realmccoi.csv'
mccoil_50_mip=[line.strip().split(',') for line in open(mccoil_50_mip)][1:]
mccoil_50_mip={'-'.join(line[0].split('-')[:2]):float(line[1]) for line in mccoil_50_mip}

mccoil_50_smarrt={sample:mccoil_100_smarrt[sample] for sample in mccoil_50_mip}

max_100_smarrt='/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_admin/seekdeep_graphs/ind_field_samples_new_replicates_new_params_389/samp_hap_counts.tsv'
max_100_smarrt=[line.strip().split('\t') for line in open(max_100_smarrt)][1:]
max_100_smarrt={line[0]:max(list(map(int, (line[1:])))) for line in max_100_smarrt}

max_50_smarrt={sample:max_100_smarrt[sample] for sample in mccoil_50_mip}

output_folder='main_figs'
subprocess.call(f'mkdir {output_folder}', shell=True)

def dict_to_lists(input_dict, type_name):
	sample_list, value_list, type_list=[],[],[]
	for sample in input_dict:
		sample_list.append(sample)
		value_list.append(input_dict[sample])
		type_list.append(type_name)
	return sample_list, value_list, type_list

def make_mccoil_df(mccoil_50_mip, mccoil_50_smarrt, mccoil_100_smarrt):
	import pandas as pd
	import scipy.stats as stats
	mip_samples, mip_values, mip_type=dict_to_lists(mccoil_50_mip, 'mccoil_50_mip')
	fifty_samples, fifty_values, fifty_type=dict_to_lists(mccoil_50_smarrt, 'mccoil_50_SMARRT')
	hundred_samples, hundred_values, hundred_type=dict_to_lists(mccoil_100_smarrt, 'mccoil_100_SMARRT')
	first_p=stats.ttest_ind(mip_values, fifty_values, equal_var=False)[1]
	second_p=stats.ttest_ind(mip_values, hundred_values, equal_var=False)[1]
	third_p=stats.ttest_ind(fifty_values, hundred_values, equal_var=False)[1]
	p_values=[first_p, second_p, third_p]
	overall_samples=mip_samples+fifty_samples+hundred_samples
	overall_values=mip_values+fifty_values+hundred_values
	overall_types=mip_type+fifty_type+hundred_type
	df=pd.DataFrame({'types': overall_types, 'values':overall_values, 'samples': overall_samples})
	return df, p_values

def make_max_df(max_50_smarrt, max_100_smarrt):
	import pandas as pd
	import scipy.stats as stats
	fifty_samples, fifty_values, fifty_type=dict_to_lists(max_50_smarrt, 'max_50_SMARRT')
	hundred_samples, hundred_values, hundred_type=dict_to_lists(max_100_smarrt, 'max_100_SMARRT')
	p_values=[stats.ttest_ind(fifty_values, hundred_values, equal_var=False)[1]]
	overall_samples=fifty_samples+hundred_samples
	overall_values=fifty_values+hundred_values
	overall_types=fifty_type+hundred_type
	df=pd.DataFrame({'types': overall_types, 'values':overall_values, 'samples': overall_samples})
	return df, p_values

def plot_it(df, title_text):
	import plotly.express as px
	fig=px.strip(df, x='types', y='values', title=title_text, hover_data=['samples']).update_traces(jitter = 1).update_traces(marker=dict(color='black'))
	dm = df.groupby('types').mean(numeric_only=True)
	ds = df.groupby('types').std(numeric_only=True)

	#use plotly express to create a scatterplot of error bars, and update all
	#traces to be red
	fig.add_scatter(x=dm.index, y=dm['values'], 
		error_y_array=ds['values'],
			mode='markers', showlegend=False, fillcolor='red')
	#fig.update_traces(marker=dict(color='red'))
#	fig.update_xaxes(title_text="")
#	fig.update_yaxes(title_text="Dd2 minor fraction (percent)")
	return fig

def add_p_value_annotation(fig, array_columns, p_values, subplot=None, _format=dict(interline=0.07, text_height=1.07, color='black')):
	''' 
	Taken from stackOverflow, here:
	https://stackoverflow.com/questions/67505252/plotly-box-p-value-significant-annotation
	Adds notations giving the p-value between two box plot data (t-test two-sided comparison)
	
	Parameters:
	----------
	fig: figure
		plotly boxplot figure
	array_columns: np.array
		array of which columns to compare 
		e.g.: [[0,1], [1,2]] compares column 0 with 1 and 1 with 2
	subplot: None or int
		specifies if the figures has subplots and what subplot to add the notation to
	_format: dict
		format characteristics for the lines

	Returns:
	-------
	fig: figure
		figure with the added notation
	'''
	# Specify in what y_range to plot for each pair of columns
	from scipy import stats
	import plotly.express as px
	import plotly.graph_objects as go
	import numpy as np
	y_range = np.zeros([len(array_columns), 2])
	for i in range(len(array_columns)):
		y_range[i] = [1.01+i*_format['interline'], 1.02+i*_format['interline']]

	# Get values from figure
	fig_dict = fig.to_dict()

	# Get indices if working with subplots
	if subplot:
		if subplot == 1:
			subplot_str = ''
		else:
			subplot_str =str(subplot)
		indices = [] #Change the box index to the indices of the data for that subplot
		for index, data in enumerate(fig_dict['data']):
			#print(index, data['xaxis'], 'x' + subplot_str)
			if data['xaxis'] == 'x' + subplot_str:
				indices = np.append(indices, index)
		indices = [int(i) for i in indices]
		print((indices))
	else:
		subplot_str = ''

	# Print the p-values
	for index, column_pair in enumerate(array_columns):
		if subplot:
			data_pair = [indices[column_pair[0]], indices[column_pair[1]]]
		else:
			data_pair = column_pair

		# Mare sure it is selecting the data and subplot you want
		#print('0:', fig_dict['data'][data_pair[0]]['name'], fig_dict['data'][data_pair[0]]['xaxis'])
		#print('1:', fig_dict['data'][data_pair[1]]['name'], fig_dict['data'][data_pair[1]]['xaxis'])

		# Get the p-value
		pvalue = p_values[index]
		#stats.ttest_ind(
		#	fig_dict['data'][data_pair[0]]['y'],
		#	fig_dict['data'][data_pair[1]]['y'],
		#	equal_var=False,
		#)[1]
		print('pvalue was', pvalue)
		if pvalue >= 0.05:
			symbol = 'ns'
		else:
			symbol = str(round(pvalue, 4))
			print('symbol was', symbol)
		# Vertical line
		fig.add_shape(type="line",
			xref="x"+subplot_str, yref="y"+subplot_str+" domain",
			x0=column_pair[0], y0=y_range[index][0], 
			x1=column_pair[0], y1=y_range[index][1],
			line=dict(color=_format['color'], width=2,)
		)
		# Horizontal line
		fig.add_shape(type="line",
			xref="x"+subplot_str, yref="y"+subplot_str+" domain",
			x0=column_pair[0], y0=y_range[index][1], 
			x1=column_pair[1], y1=y_range[index][1],
			line=dict(color=_format['color'], width=2,)
		)
		# Vertical line
		fig.add_shape(type="line",
			xref="x"+subplot_str, yref="y"+subplot_str+" domain",
			x0=column_pair[1], y0=y_range[index][0], 
			x1=column_pair[1], y1=y_range[index][1],
			line=dict(color=_format['color'], width=2,)
		)
		## add text at the correct x, y coordinates
		## for bars, there is a direct mapping from the bar number to 0, 1, 2...
		fig.add_annotation(dict(font=dict(color=_format['color'],size=14),
			x=(column_pair[0] + column_pair[1])/2,
			y=y_range[index][1]*_format['text_height'],
			showarrow=False,
			text=symbol,
			textangle=0,
			xref="x"+subplot_str,
			yref="y"+subplot_str+" domain"
		))
		fig.update_layout(margin=dict(l=30, r=20, t=90, b=20))
	return fig

mccoil_df, mccoil_p_values=make_mccoil_df(mccoil_50_mip, mccoil_50_smarrt, mccoil_100_smarrt)
mccoil_df.to_csv(f'{output_folder}/fig2_field_mccoil_COI.tsv', sep='\t', index=False)
mccoil_fig=plot_it(mccoil_df, 'Real McCoil COI')
mccoil_fig=add_p_value_annotation(mccoil_fig, [[0,1],[0,2],[1,2]], mccoil_p_values)
mccoil_fig.write_html(f'{output_folder}/fig2_field_mccoil_COI.html')
mccoil_fig.write_image(f'{output_folder}/fig2_field_mccoil_COI.svg')
mccoil_fig.write_image(f'{output_folder}/fig2_field_mccoil_COI.png')


max_df, max_p_values=make_max_df(max_50_smarrt, max_100_smarrt)
max_df.to_csv(f'{output_folder}/fig2_field_max_COI.tsv', sep='\t', index=False)
max_fig=plot_it(max_df, 'Max # Haplotypes')
max_fig=add_p_value_annotation(max_fig, [[0,1]], max_p_values)
max_fig.write_html(f'{output_folder}/fig2_field_max_COI.html')
max_fig.write_image(f'{output_folder}/fig2_field_max_COI.svg')
max_fig.write_image(f'{output_folder}/fig2_field_max_COI.png')
