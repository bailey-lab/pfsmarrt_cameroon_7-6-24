'''
creates the 'population haplotypes' and 'heterozygosity' columns of table 4
'''
import gzip
prefix=('/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_admin/'
	'seekdeep_outputs/field_samples_new_replicates_new_params_389/seekdeep_'
	'output_round2/analysis/popClustering/')
suffix='/analysis/selectedClustersInfo.tab.txt.gz'

amplicons=['ama1', 'heome-a', 'heome-b', 'heome-c', 'heome-d', 'heome-e',
'heome-f', 'heome-g', 'heome-h']

def calc_hetero(hap_dict):
	'''
	heterozygosity is the probability that two alleles chosen with probability
	equal to the population allele frequency will be different from each other.
	This is 1-(sum of probabilities of drawing the same allele twice). sum of
	probabilities of drawing the same allele twice is P(allele_1)^2+
	P(allele_2)^2...+P(allele_n)^2 where n is total number of alleles (and
	allele count is same as haplotype count)
	'''
	from decimal import Decimal
	homozygous_prob=Decimal('0.0')
	for hap in hap_dict:
		homozygous_prob+=Decimal(hap_dict[hap])**2
	return round(float(Decimal('1.0')-homozygous_prob), 3)

for amplicon in amplicons:
	path=prefix+amplicon+suffix
	line_count=0
	hap_dict={}
	for line_number, line in enumerate(gzip.open(path, mode='rt')):
		line=line.strip().split('\t')
		if line_number==0:
			header_dict={}
			for column_number, column in enumerate(line):
				header_dict[column]=column_number
		else:
			hap=line[header_dict['h_popUID']]
			frac=line[header_dict['h_PopFrac']]
			if hap not in hap_dict:
				hap_dict[hap]=frac
			elif frac!=hap_dict[hap]:
				print('in line number', line_number, 'getting', hap, frac, 'but expected', hap_dict[hap])
	heterozygosity=calc_hetero(hap_dict)
	print(amplicon, len(hap_dict), heterozygosity)