import yaml
mapped_haplotypes=snakemake.input.mapped_haplotypes
count_dict=yaml.safe_load(open(snakemake.input.counts))
parsed_haplotypes=open(snakemake.output.parsed_haplotypes, 'w')
all_haplotypes=open(snakemake.output.all_haplotypes, 'w')

def make_hap_dict(mapped_haplotypes):
	'''
	gets all the genomes that best associated with each haplotype (as well as
	their scores and the lengths of each haplotype). format is:
	{hap:[score, {genomes}, hap_size]}
	'''
	hap_dict={}
	for haplotype_file in mapped_haplotypes:
		genome=haplotype_file.split('/')[-1].replace('_blatted_haplotypes.psl', '')
		print(haplotype_file)
		for line in open(haplotype_file):
			line=line.strip().split('\t')
			if line[0].isdigit() and len(line)>3:
				score, hap, hap_size=line[0], line[9], line[10]
				hap=hap.split('_')[0]
				if hap not in hap_dict or int(score)>hap_dict[hap][0]:
					hap_dict[hap]=[int(score), {genome}, int(hap_size)]
				elif int(score)==hap_dict[hap][0]:
					hap_dict[hap][1].add(genome)
	return hap_dict

hap_dict=make_hap_dict(mapped_haplotypes)

print('hap dict is', hap_dict)

printing_dict={}
for amplicon in count_dict:
	for sample in count_dict[amplicon]:
		for replicate in count_dict[amplicon][sample]:
			printing_dict.setdefault(replicate, {})
			printing_dict[replicate].setdefault(amplicon, {})
			total, hap_count=0,0
			for hap in count_dict[amplicon][sample][replicate]:
				count=count_dict[amplicon][sample][replicate][hap]
				total+=count
				hap_count+=1
			printing_dict[replicate][amplicon]['hap_count']=hap_count
			for hap in count_dict[amplicon][sample][replicate]:
				count=count_dict[amplicon][sample][replicate][hap]
				score, genomes, hap_size=hap_dict[hap]
				for genome in genomes:
					printing_dict[replicate][amplicon].setdefault(genome, [])
					printing_dict[replicate][amplicon][genome].append([hap, score, hap_size, count, total])
parsed_haplotypes.write('replicate\tamplicon\thap_count\tPf3D7\tPfDd2\tPf7G8\n')
all_haplotypes.write('replicate\tamplicon\thap_count\tPf3D7\tPfDd2\tPf7G8\n')
all_genomes=['Pf3D7', 'PfDd2', 'Pf7G8']
for replicate in printing_dict:
	printing_line=[replicate]
	for amplicon in printing_dict[replicate]:
		genomes=set(printing_dict[replicate][amplicon].keys())-{'hap_count'}
		hap_count=printing_dict[replicate][amplicon]['hap_count']
		bad_rep=False
		printing_line=[replicate, amplicon, str(hap_count)]
		for genome in genomes:
			printing_list=[]
			for hap_list in printing_dict[replicate][amplicon][genome]:
				hap, score, hap_size, count, total=hap_list
				string_version=f'{hap} gen_match: {score}/{hap_size} replicate_frac: {count}/{total} ({round(count/total*100, 1)}%)'
				printing_list.append(string_version)
				if hap_count>1 and '955' not in replicate:
					bad_rep=True
				elif hap_count==1 and 'DD2' in replicate and 'PfDd2' not in genomes:
					bad_rep=True
				elif hap_count==1 and 'DD2' in replicate and genome=='PfDd2' and score/hap_size<1:
					bad_rep=True
				elif hap_count==1 and '7G8' in replicate and 'Pf7G8' not in genomes:
					bad_rep=True
				elif hap_count==1 and '7G8' in replicate and genome=='Pf7G8' and score/hap_size<1:
					bad_rep=True
				elif '7G8' not in replicate and 'DD2' not in replicate and hap_count==1 and 'Pf3D7' not in genomes:
					bad_rep=True
				elif hap_count==1 and genome=='Pf3D7' and score/hap_size<1:
					bad_rep=True
				elif hap_count==1 and '955' in replicate:
					bad_rep=True
				elif hap_count>2 and '955' in replicate:
					bad_rep=True
				elif hap_count>1 and '955' in replicate and genome=='Pf3D7' and count/total<0.85:
					bad_rep=True
				elif hap_count>1 and '955' in replicate and genome=='Pf3D7' and score/hap_size<1:
					bad_rep=True
			printing_dict[replicate][amplicon][genome]=printing_list
		for genome in all_genomes:
			if genome in genomes:
				printing_list=printing_dict[replicate][amplicon][genome]
				printing_line.append(' *** '.join(printing_list))
			else:
				printing_line.append('no match')
		if bad_rep:
			parsed_haplotypes.write('\t'.join(printing_line)+'\n')
		all_haplotypes.write('\t'.join(printing_line)+'\n')