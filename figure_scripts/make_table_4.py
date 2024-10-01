'''
creates the 'population haplotypes' column of table 4
'''
import gzip
prefix=('/nfs/jbailey5/baileyweb/asimkin/pf_SMARRT_methods/alfred_admin/'
	'seekdeep_outputs/field_samples_new_replicates_new_params_389/seekdeep_'
	'output_round2/analysis/popClustering/')
suffix='/analysis/population/PopSeqs.fastq.gz'

amplicons=['ama1', 'heome-a', 'heome-b', 'heome-c', 'heome-d', 'heome-e',
'heome-f', 'heome-g', 'heome-h']

for amplicon in amplicons:
	path=prefix+amplicon+suffix
	line_count=0
	for line in gzip.open(path, mode='rt'):
		if len(line)>1:
#			print('line is', line)
			line_count+=1
	print(amplicon, line_count/4)