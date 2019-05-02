import re

def fasta_reader(filename):
	'''
	Read in fasta file as dictionary, with headers as keys and 
	sequences as values
	'''
	infile = open(filename)

	seqs = {}

	while 1:
		line = infile.readline()
		if not line:
			break
		# skip first character '>'
		header = line.strip()[1:]
		seq = infile.readline().strip()
		seqs[header] = seq

	return seqs

def fasta_writer(outfile, seqs):
	for x in seqs:
		outfile.write('>'+x+'\n'+seqs[x]+'\n')
	outfile.close()

def check_REs(lib):
	'''
	Remove any sequences that contain restriction sites used for cloning
	'''
	good_seqs = {}

	# just check by hand, easier
	# AleI, CACNNNNGTG
	# MluI, ACGCGT
	# KpnI, GGTACC
	# XbaI, TCTAGA
	# SpeI, ACTAGT

	patterns = [re.compile('CAC[ACTG]{4}GTG'), re.compile('ACGCGT'),
				re.compile('GGTACC'), re.compile('TCTAGA'), re.compile('ACTAGT')]


	for x in lib:
		# trim sequences so as not to include the KpnI and MluI site that were added in
		matches = [re.search(pattern, lib[x][25:-25]) for pattern in patterns]
		match_count = sum([1 for match in matches if match != None])
		if match_count == 0:
			good_seqs[x] = lib[x]
		else:
			print(x, "removed")

	print(len(lib) - len(good_seqs), "sequences removed.")
	print(len(good_seqs), "sequences after removing those with restriction sites.")
	return good_seqs
	

def reverse_complement(seq):
	"""
	Return the reverse complement of a nucleotide string
	"""
	complement = {'A': 'T', 'T':'A', 'C':'G', 'G':'C'}
	
	rc = ''.join([complement[nt] for nt in seq[::-1]])
	return rc


def best_A_content(oligo):
	'''
	Choose the strand with the lowest A content because A's are harder to
	synthesize.
	'''

	# get the reverse compliment of the oligo
	rc_oligo = reverse_complement(oligo)

	# count the number of A's for both strands
	oligo_As = sum( [1 for nt in oligo if nt == 'A'] )
	rc_As = sum( [1 for nt in rc_oligo if nt == 'A'] )

	# save the oligo with lower number of A's
	if oligo_As < rc_As:
		final_oligo = oligo
	else:
		final_oligo = rc_oligo

	return final_oligo

def csv_writer(lib, filename):
	with open(filename, 'w') as outfile:
		for x in lib:
			outfile.write(x+','+lib[x]+'\n')


if __name__ == '__main__':

	# read in primers
	fwd_primers = fasta_reader('fwd_primers_cre.fasta')
	# reverse primer assumed to be given in 3' to 5' on bottom strand,
	# need to reverse complement
	rev_primers = fasta_reader('rev_primers_cre.fasta')
	rev_primers = {x : reverse_complement(rev_primers[x]) for x in rev_primers}

	# read in enhancer templates for variant libraries
	enhancers = fasta_reader('enhancer_templates.fasta')

	# dictionary that will hold dictionaries of the subpool libraries
	# key will be name of subpool, easy to match with primers
	subpool_libraries = {}
#-------------------------------------------------------------------------------
	# subpool 2, distance to promoter
	dist_to_promoter = {}
	# putting these in lowercase makes it easier to double check work
	# will convert to all uppercase for final version to be sent to Agilent
	consensus = 'TGACGTCA'.lower()
	consensus_flank = 'ATTGACGTCAGC'.lower()

	for x in enhancers:
		enh = enhancers[x]
		seqs = [ enh[:i] + consensus + enh[i+len(consensus):] 
					for i in range(len(enh)-len(consensus)+1)]

		# add to dictionary with appropriate headers
		consensus_seqs = {'consensus_dist_' + str(142-i) + '_' + x :  seqs[i] 
								for i in range(len(seqs))}

		dist_to_promoter.update(consensus_seqs)

		# now do same thing for consensus+flank
		seqs = [ enh[:i] + consensus_flank + enh[i+len(consensus_flank):] 
					for i in range(len(enh)-len(consensus_flank)+1)]

		# add to dictionary with appropriate headers
		consensus_flank_seqs = {'consensusflank_dist_' + str(138-i) + '_' + x :  seqs[i] 
									for i in range(len(seqs))}
		
		dist_to_promoter.update(consensus_flank_seqs)

	subpool_libraries['subpool_2'] = dist_to_promoter
	print(len(dist_to_promoter), "sequences in subpool 2")

#-------------------------------------------------------------------------------
	var_cre_sites = {}

	# subpool 4, variable distance between two cre sites using combinations of 
	# consensus, weak, moderate and half flank sites. 
	weak = 'ATTGAAGTCAGC'.lower()
	moderate = 'ATTGACGTCTGC'.lower()
	half_flank = 'XXXGCCGTCATA'.lower()
	cre_vars = [consensus_flank, weak, moderate, half_flank]
	# used for headers
	names = ['consensus', 'weak', 'moderate', 'half']
	spacing = [0, 1, 6, 11, 16]


	for x in enhancers:
		enh = enhancers[x]
		for j in range(len(cre_vars)):
			for k in range(len(cre_vars)):
				if names[j] == names[k] == 'half': # don't include half X half site to save space
					continue

				site1 = cre_vars[j]
				site2 = cre_vars[k]

				for space in spacing:
					# 0 bp spacing, flank + site1 + site2 + flank
					# trim last two bases from site1 and first two bases from site 2
					if space == 0:
						seqs = [ enh[:i] + site1[:-2] + site2[2:] + enh[i+20:]
							for i in range(0, len(enh)-20 + 1, 5)]
						var_cre = {names[j]+'_'+names[k]+'0 bp spacing_dist_'+str((150 - 20-(i*5))) + '_' + x : seqs[i]
									for i in range(len(seqs))}
					else:
						seqs = [ enh[:i] + site1 + enh[i+len(site1) : i+len(site1)+space] + site2 + enh[i+2*len(site2)+space:]
							for i in range(0, len(enh)-(2*len(site1) + space) + 1, 5)]

						var_cre = {names[j]+'_'+names[k]+ '_' + str(space)+ 'bp spacing_dist_'+str((150 - (2*len(consensus_flank)+space)-(i*5))) + '_' + x : seqs[i]
										for i in range(len(seqs))}

					# go through sites and find length of string of x's and
					# replace with appropriate bases from template
					pattern = re.compile('x+') # one or more x's
					for key in var_cre:
						oligo = var_cre[key]
						loc = oligo.find('x')
						if loc > -1 :
							num_x = len(max(pattern.findall(oligo)))
							# this is the half site on the right, instead of having
							# one base in between being the template assign it to
							# G from the other half site
							if num_x == 1:
								new_oligo = oligo.replace('x', 'G')
							else:
								new_oligo = oligo[:loc] + enh[loc:loc+num_x] + oligo[loc+num_x:]
							var_cre[key] = new_oligo
					
					var_cre_sites.update(var_cre)

	subpool_libraries['subpool_4'] = var_cre_sites
	print(len(var_cre_sites), "sequences in subpool 4")

#-------------------------------------------------------------------------------

	# subpool 5, 6 combinations of no addition, consensus and weak site, all with 
	# flanks
	cre_vars = ['xxxxxxxxxxxx', consensus_flank, weak]
	names = ['no_site', 'consensus', 'weak']
	cre_combs = {}

	for x in enhancers:
		enh = enhancers[x]
		# the six sites are in a constant position on the template, so parse out
		# the parts of the template we're going to use now
		spacers = [enh[i:i+13] for i in range(0, len(enh), 25)]
		for i in range(len(cre_vars)):
			for j in range(len(cre_vars)):
				for k in range(len(cre_vars)):
					for l in range(len(cre_vars)):
						for m in range(len(cre_vars)):
							for n in range(len(cre_vars)):
								seq = spacers[0] + cre_vars[i] + spacers[1] + cre_vars[j]
								seq += spacers[2] + cre_vars[k] + spacers[3] + cre_vars[l]
								seq += spacers[4] + cre_vars[m] + spacers[5] + cre_vars[n] 

								# replace all  xxx no site with bases from template
								loc = seq.find('xxxxxxxxxxxx')
								while loc > -1:
									seq = seq[:loc] + enh[loc:loc+12] + seq[loc+12:]
									loc = seq.find('xxxxxxxxxxxx')

								header = '_'.join([names[i], names[j], names[k], names[l], names[m], names[n], x])
								cre_combs[header] = seq
	subpool_libraries['subpool_5'] = cre_combs
	print(len(cre_combs), "sequences in subpool 5")
    
    #-------------------------------------------------------------------------------
	# subpool 3, variable distance between two consensus-flank sequences
	consensus_flank_spacing = {}

	# 0 bp spacing
	flank_con_con_flank = 'ATTGACGTCATGACGTCAGC'.lower()
	for x in enhancers:
		enh = enhancers[x]
		seqs = [ enh[:i] + flank_con_con_flank + enh[i+len(flank_con_con_flank):] 
					for i in range(len(enh)-len(flank_con_con_flank)+1)]
		twobs_seqs = {'2BS 0 bp spacing flank-con-con-flank_dist_' + str(130-i) + '_' + x : seqs[i]
							for i in range(len(seqs))}
		consensus_flank_spacing.update(twobs_seqs)

	# this refers to number of bases in between the two sites
	spacing = [1, 6, 11, 16, 66]
	for x in enhancers:
		enh = enhancers[x]
		for space in spacing:
			seqs = [ enh[:i] + consensus_flank + enh[i+len(consensus_flank) : i+len(consensus_flank)+space] + consensus_flank + enh[i+2*len(consensus_flank)+space:]
						for i in range(len(enh)-(2*len(consensus_flank) + space) + 1)]

			twobs_spacing = {'2BS ' + str(space) + ' bp spacing consensus+flank x2_dist_' + str((150 - (2*len(consensus_flank)+space)-i)) + '_' + x : seqs[i]
								for i in range(len(seqs))}
			consensus_flank_spacing.update(twobs_spacing)

	subpool_libraries['subpool_3'] = consensus_flank_spacing
	print(len(consensus_flank_spacing), "sequences in subpool 3")

	#-------------------------------------------------------------------------------
	# output rough draft of lib to check by eye, templates in upper case and
	# sites in lower case
	full_lib = {}
	for x in subpool_libraries:
		full_lib.update(subpool_libraries[x])

	csv_writer(full_lib, 'cre_lib_rough_draft.csv')

	# read in controls
	controls = fasta_reader('cre_controls.fasta')

	full_lib = {}

	MluI = 'ACGCGT'
	KpnI = 'GGTACC'

	# add controls to each subpool, then add corresponding primers
	for x in subpool_libraries:
		subpool = subpool_libraries[x]
		subpool.update(controls)
		fwd_primer = fwd_primers[x+'_F']
		rev_primer = rev_primers[x+'_R']
		for name in subpool:
			full_seq = fwd_primer + MluI + subpool[name] + KpnI + rev_primer
			subpool[name] = full_seq.upper()
		full_lib.update(subpool)

	print(len(full_lib), "sequences in library before checking REs")
	cleaned_lib = check_REs(full_lib)

	print("Selecting strand with lower A content")
	lib_final = {x : best_A_content(cleaned_lib[x]) for x in cleaned_lib}

	csv_writer(lib_final, 'cre_lib_final.csv')













