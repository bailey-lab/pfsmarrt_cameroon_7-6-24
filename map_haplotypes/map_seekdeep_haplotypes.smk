configfile: 'map_seekdeep_haplotypes.yaml'

rule all:
	input:
		parsed_haplotypes=config['output_folder']+'/unexpected_haplotypes.tsv'

rule gather_haplotypes:
	'''
	gathers all haplotype sequences in a single fastq file
	'''
	input:
		analysis_folder=config['original_analysis_folder']
	output:
		haplotypes=config['output_folder']+'/all_haplotypes.fasta'
	script:
		'scripts/gather_haplotypes.py'

rule map_haplotypes:
	'''
	maps haplotypes to genomes using blat
	'''
	input:
		haplotypes=config['output_folder']+'/all_haplotypes.fasta',
		genome=config['genomes']+'/{genome}.fasta'
	output:
		mapped_haplotypes=config['output_folder']+'/{genome}_blatted_haplotypes.psl'
	shell:
		'blat {input.genome} {input.haplotypes} {output.mapped_haplotypes}'

rule parse_blat:
	'''
	parses blat output files to see which genomes each haplotype maps best to
	'''
	input:
		mapped_haplotypes=expand(config['output_folder']+'/{genome}_blatted_haplotypes.psl', genome=['Pf3D7', 'PfDd2', 'Pf7G8']),
		counts=config['count_file']
	output:
		all_haplotypes=config['output_folder']+'/all_haplotypes.tsv',
		parsed_haplotypes=config['output_folder']+'/unexpected_haplotypes.tsv'
	script:
		'scripts/parse_blat.py'