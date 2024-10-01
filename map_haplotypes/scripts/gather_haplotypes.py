import os
import gzip
analysis_folder=snakemake.input.analysis_folder
haplotypes=open(snakemake.output.haplotypes, 'w')
clustering_folder=analysis_folder+'/popClustering'

folders=set(os.listdir(clustering_folder))-{'locationByIndex'}
for folder in folders:
	fastq_file=f'{clustering_folder}/{folder}/analysis/population/PopSeqs.fastq.gz'
	for line_number, line in enumerate(gzip.open(fastq_file, mode='rt')):
		if line_number%4==0:
			haplotypes.write('>'+line[1:])
		elif line_number%4==1:
			haplotypes.write(line)