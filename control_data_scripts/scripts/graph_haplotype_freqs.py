import yaml
import plotly.express as px
import pandas as pd

yaml_file='/nfs/jbailey5/baileyweb/asimkin/other_people/jmsadler/bigv10ctrl/analyzed_by_PCR_replicates/seekdeep_analysis_snakemake/half_percent_exact_only_analysis_graphs/read_dict.yaml'
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
				new_sample=sample.replace('-', '.')
				hap_dict={}
				for replicate in main_db[amplicon][sample]:
					for hap in main_db[amplicon][sample][replicate]:
						hap_count=main_db[amplicon][sample][replicate][hap]
						hap_dict[hap]=hap_dict.setdefault(hap, 0)+hap_count
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

hap_counts=reorganize_db(main_db)
frac_counts=get_fracs(hap_counts)
graphing_dict=flatten_graphing(frac_counts)
for amplicon in graphing_dict:
	output_path='stacked_bars/'+amplicon+'_stacked_bar.html'
	df=graphing_dict[amplicon]
	fig=px.bar(df, x="samples", y="fracs", color="haps", title=amplicon)
	fig.write_html(output_path)
