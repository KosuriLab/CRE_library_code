import re
import random
import numpy as np


random.seed(123)


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

def csv_reader(filename):
	infile = open(filename)
	seqs = {}

	for line in infile:
		name, seq = line.strip().split(',')
		seqs[name] = seq

	return seqs

def fasta_writer(filename, seqs):
	outfile = open(filename, 'w')
	
	for x in seqs:
		outfile.write('>'+x+'\n'+seqs[x]+'\n')

	outfile.close()

def check_REs(lib):
	'''
	Remove any sequences that contain restriction sites used for cloning
	'''
	good_seqs = {}

	# just check by hand, easier
	# MluI, ACGCGT
	# KpnI, GGTACC
	# XbaI, TCTAGA
	# SpeI, ACTAGT

	patterns = [re.compile('ACGCGT'), re.compile('GGTACC'), 
		re.compile('TCTAGA'), re.compile('ACTAGT')]


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
            
def capitalize_values(d):
    result = {}
    for key, value in d.items():
        upper_value = value.upper()
        result[key] = upper_value
    return result

def similar(a, b):
    # must be same length sequence
    num_same = sum([1 for i in range(len(a)) if a[i] == b[i]])
    return num_same

def most_scramble(sequence, n_iter):
    '''
    Scramble sequence n_iter times and choose the most scrambled sequence
    '''
    scrambled_seqs = [''.join(random.sample(list(sequence), 
                                            len(sequence))) for i in range(n_iter)]
    similarity = [similar(sequence, x) for x in scrambled_seqs]
    max_scrambled = scrambled_seqs[np.argmin(similarity)]
    return max_scrambled

def calculate_gc(seq):
	count = 0
	for i in range(len(seq)):
		if seq[i] == 'G' or seq[i] == 'C' or seq[i] == 'g' or seq[i] == 'c':
			count += 1

	return count / float(len(seq))

def name_background(backgrounds):
    bg_names = []
    name_number = 1
    seq = list(backgrounds.keys())
    seq_name = {}
    
    for bg in backgrounds.values():
        bg_name = 'background_' + str(name_number)
        name_number = name_number + 1
        bg_names.append(bg_name)
        
    for i in range(len(seq)):
        seq_name[bg_names[i]] = seq[i]
        
    return(seq_name)


#------------------------------------------------------------------------
if __name__ == '__main__':
    
    #import old backgrounds to be used in some designs
    old_bgs = csv_reader('oldbackgrounds.csv')
    
    #Generating the gc content background library with CREs
    bgs_top3 = capitalize_values(csv_reader('bgs_top3.csv'))
    
    #import 3 positive CRE controls
    controls = csv_reader('controls.csv')
    
    full_lib = {}

#------------------------------------------------------------------------------
#Add name for backgrounds with number
    
bgs_top3_name = name_background(bgs_top3)

combined_bgs = {**bgs_top3_name, **old_bgs}


csv_writer(bgs_top3_name, 'bgs_top3_name.csv')


#Scramble 10 bp of the backgrounds in 5 bp increments to check for TFBSs 
#generated

tiles = {}

for x in combined_bgs:
    
    bg = combined_bgs[x]
    counter = 0
    re_patterns = [re.compile('ACGCGT'), re.compile('GGTACC'), 
                   re.compile('TCTAGA'), re.compile('ACTAGT'),
                   re.compile('CGTCA'), re.compile('TGACG')]
    
    for i in range(0, len(bg), 5):
        if i + 10 > len(bg):
            continue
        to_scramble = bg[i:i+10]
        scrambled = most_scramble(to_scramble, 100)
        tile = bg[:i] + ''.join(scrambled).lower() + bg[i+10:]
        matches = [re.search(pattern, tile.upper()) for pattern in re_patterns]
        match_count = sum([1 for match in matches if match != None])
        if match_count == 0:
            kept_scramble = tile
        else:
            scrambled = most_scramble(to_scramble, 1000)
            tile = bg[:i] + ''.join(scrambled).lower() + bg[i+10:]
            matches = [re.search(pattern, tile.upper()) for pattern in re_patterns]
            match_count = sum([1 for match in matches if match != None])
            if match_count == 0:
                kept_scramble = tile
            else:
                scrambled = most_scramble(to_scramble, 10)
                tile = bg[:i] + ''.join(scrambled).lower() + bg[i+10:]
                matches = [re.search(pattern, tile.upper()) for pattern in re_patterns]
                match_count = sum([1 for match in matches if match != None])
                if match_count == 0:
                    kept_scramble = tile
                else:
                    counter = 1 + counter
                    
        tile_name = 'scramble_dist_' + str(140-i) + '_similarity_' + str(similar(to_scramble, scrambled)) + '_' + x
        tiles[tile_name] = kept_scramble
            
    full_lib.update(tiles)
    print(len(tiles), 'scrambled sequences in minP subpool')
    print(str(counter), ' scrambled sequences still contain unwanted motifs')
    
    
csv_writer(tiles, 'cre_test_scramble.csv')
           

#------------------------------------------------------------------------------
#Add 0-6 consensus CRE sites to backgrounds, use old backgrounds as a control

#putting these in lowercase makes it easier to double check work will convert 
#to all uppercase for final version to be sent to Agilent

consensus = 'TGACGTCA'.lower()
consensus_flank = 'ATTGACGTCAGC'.lower()
nosite = 'XXXXXXXX'.lower()
nosite_flank = 'XXXXXXXXXXXX'.lower()

cre_vars = [nosite_flank, consensus_flank]
names = ['nosite', 'consensus']
cre_combs = {}
bg = []

for x in combined_bgs:
    bg = combined_bgs[x]
    # the six sites are in a constant position on the template, so parse out
    # the parts of the template we're going to use now
    spacers = [bg[i:i+13] for i in range(0, len(bg), 25)]
    
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
                                seq = seq[:loc] + bg[loc:loc+12] + seq[loc+12:]
                                loc = seq.find('xxxxxxxxxxxx')
                                
                            gc = calculate_gc(seq)
                            
                            header = 'sixsite_' + '_'.join([names[i], names[j], 
                                       names[k], names[l], names[m], names[n], 
                                       x, 'gc', str(gc)])
        
                            cre_combs[header] = seq
                            full_lib.update(cre_combs)
                            
        print(len(cre_combs), "sixsite sequences in minP subpool")
        

csv_writer(cre_combs, 'cre_test_sixsite.csv')

#------------------------------------------------------------------------------
#Design the 2-site library again with CREs spaced 5-20 bp apart and CRE sites 
#varying between consensus CREs and their reverse complement. Move both CREs
#along backgrounds maintaining spacing. Keep in mind both CRE sites are 12 bp
#long with flanks.

cre_space_dist = {}
bg = []
seqs = []
spacing = []
lessthan5bp_seqs = {}
seq_name = {}

#generate special cases where full flanks aren't used for spacing < 5 bp

lessthan5bpspacing = ['ATTGACGTCATTGACGTCAGC'.lower(),
                      'ATTGACGTCAGTTGACGTCAGC'.lower(),
                      'ATTGACGTCAGATTGACGTCAGC'.lower(),
                      'ATTGACGTCAGCATTGACGTCAGC'.lower()]
lessthan5bpnames = ['space_1_', 'space_2_', 'space_3_', 'space_4_']

for x in old_bgs:
    bg = old_bgs[x]
    
    for j in range(len(lessthan5bpspacing)):
        spacing = lessthan5bpspacing[j]
        
        seqs = [bg[:i] + spacing + bg[i+len(spacing):] for i in range(len(bg)-len(spacing)+1)]
        
        #distance is written accommodating the actual 2 bp + distance to the proximal CRE
        
        lessthan5bp_seqs = {'twosite_' + lessthan5bpnames[j] + 'dist_' + 
                            str(150 + 2 - len(spacing) - i) + '_' + 
                            x : seqs[i] for i in range(len(seqs))}
        
        cre_space_dist.update(lessthan5bp_seqs)
        
    print(len(cre_space_dist), "twosite sequences in minP subpool")
        
#generate spacing starting at a total of 5 and ending at 13 to accomodate a
#full translocation of CRE from 13 to 21 bp along the DNA. Spacing below refers 
#to number of bases in between the two CREs with flanks so spacing + 4 refers 
#to the actual number of bases between CREs

for x in old_bgs:
    bg = old_bgs[x]
    for space in range(1, 10, 1):
        seqs = [bg[:i] + consensus_flank + bg[i + 12 : i + 12 + space] + consensus_flank + bg[i + 2*12 + space:] for i in range(len(bg)-(2*12 + space) + 1)]
        
        #Spacing is written accommodating the actual 4 bp + spacing here between CREs
        
        seq_name = {'twosite_space_' + str(space + 4) + '_dist_' + 
                    str((150 - (2*12 + space) + 2 - i)) + '_' +
                    x : seqs[i] for i in range(len(seqs))}
        
        cre_space_dist.update(seq_name)
                
    full_lib.update(cre_space_dist)
    print(len(cre_space_dist), "twosite sequences in minP subpool")
    
    
    csv_writer(cre_space_dist, 'cre_test_twosite.csv')
  
    
#------------------------------------------------------------------------------
#2-site library design with fixed proximal CREs and moving distal CREs.
#Proximal and distal CREs will be offset by 0, 5, 10, and 15 proximal to the 
#promoter while the distal CRE moves from 0 bp spacing to the end of the
#background
    
twosite_offset_distdistal = {}
twosite_offset_distproximal = {}
bg = []
seqs = []
spacing = []
lessthan5bp_distal = {}

offsets = [0, 5, 10, 15]

#make 1-4 bp spacing separate from the rest of the library
for x in old_bgs:
    bg = old_bgs[x]
    
    for j in range(len(lessthan5bpspacing)):
        spacing = lessthan5bpspacing[j]
        
        for offset in offsets: 
            
            seqs = [bg[:len(bg) - len(spacing) - offset] + spacing + bg[len(bg) - offset:]]
            
            #distance is written accommodating the actual 2 bp + distance to the proximal CRE
            
            lessthan5bp_distal = {'twositedistdistal_offset_' + str(offset) + 
                                   '_' + lessthan5bpnames[j] + 'dist_' + 
                                   str(offset + len(spacing) - 10) + '_' + 
                                   x : seqs[i] for i in range(len(seqs))}
            
            twosite_offset_distdistal.update(lessthan5bp_distal)
            
    print(len(twosite_offset_distdistal), "twosite distal distance sequences")
            
#make all other spacings from 5 to 31 with offsets
distal = {}

for x in old_bgs:
    bg = old_bgs[x]
    
    for offset in offsets:
        seqs = [bg[:80 + i] + consensus_flank + bg[80 + i + 12:len(bg) - 12 - offset] + consensus_flank + bg[len(bg) - offset:] 
                for i in range(70 - offset - 24)]
                
        distal = {'twositedistdistal_offset_' + str(offset) + '_space_' +
                  str(70 - offset - 2*10 - i) + '_dist_' + 
                  str(150 - 80 - 12 - i + 2) + '_' +
                  x : seqs[i] for i in range(len(seqs))}
        
        twosite_offset_distdistal.update(distal)
        
    full_lib.update(twosite_offset_distdistal)
    print(len(twosite_offset_distdistal), "twosite distal distance sequences")
            
    csv_writer(twosite_offset_distdistal, 'cre_test_twosite_offset_distdistal.csv')


#------------------------------------------------------------------------------
#2-site library design with fixed distal CREs and moving proximal CREs.
#Proximal and distal CREs will be offset by 0, 5, 10, and 15 from the 5' end of
#the background while the proximal CRE moves from 0 bp spacing to the 3' end of
#the background
    
twosite_offset_distproximal = {}
bg = []
seqs = []
spacing = []
lessthan5bp_distproximal = {}

#make 1-4 bp spacing separate from the rest of the library
for x in old_bgs:
    bg = old_bgs[x]
    
    for j in range(len(lessthan5bpspacing)):
        spacing = lessthan5bpspacing[j]
        
        for offset in offsets: 
            
            seqs = [bg[:80 + offset] + spacing + bg[80 + offset + len(spacing):]]
            
            #distance is written accommodating the actual 2 bp + distance to the proximal CRE
            lessthan5bp_proximal = {'twositedistproximal_offset_' + str(offset) + 
                                   '_' + lessthan5bpnames[j] + 'dist_' + 
                                   str(150 - 80 - offset - len(spacing) + 2) + '_' + 
                                   x : seqs[i] for i in range(len(seqs))}
            
            twosite_offset_distproximal.update(lessthan5bp_proximal)
    
    print(len(twosite_offset_distproximal), "twosite proximal distance sequences")
            
#make all other spacings from 5 to 31 with offsets
proximal = {}

for x in old_bgs:
    bg = old_bgs[x]
    
    for offset in offsets: 
        
        seqs = [bg[:80 + offset] + consensus_flank + bg[80 + offset + 12:80 + offset + 12 + i + 1] + consensus_flank + bg[80 + offset + 24 + i + 1:]
                for i in range(70 - offset - 24)]
        
        #distance is written accommodating the actual 2 bp + distance to the proximal CRE
        proximal = {'twositedistproximal_offset_' + str(offset) + '_space_' + 
                    str(i + 1 + 4) + '_dist_' + str(69 - offset - 12 - 10 - i) + '_' +
                    x : seqs[i] for i in range(len(seqs))}
        
        twosite_offset_distproximal.update(proximal)
        
    full_lib.update(twosite_offset_distproximal)
    print(len(twosite_offset_distproximal), "twosite proximal distance sequences")
            
    csv_writer(twosite_offset_distproximal, 'cre_test_twosite_offset_distproximal.csv')
    

#------------------------------------------------------------------------------
#make combined library with controls, restriction enzyme sites and primer seqs

MluI = 'ACGCGT'
KpnI = 'GGTACC'

fwd_primer = 'ATTACCATGTTATCGGGCG'
rev_primer = 'CTGAGGTCGTGGTTTAGAT'

full_lib.update(controls)

final_lib = {}
counter = 0

for x in full_lib:
    sequence = full_lib[x]
    final_seq = fwd_primer + MluI + sequence + KpnI + rev_primer
    
    length = len(final_seq)
    
    if length != 200:
        counter = counter + 1
    
    final_lib[x] = final_seq
    
print('sequences not 200 nt long: ', str(counter))
    
csv_writer(final_lib, 'cre_test_final_lib.csv')

print(len(final_lib), "sequences in library before checking REs")
cleaned_lib = check_REs(capitalize_values(final_lib))

csv_writer(cleaned_lib, 'cre_test_cleaned_lib.csv')

print("Selecting strand with lower A content")
Acontent_lib = {x : best_A_content(cleaned_lib[x]) for x in cleaned_lib}

fasta_writer('20181220_cre_Acontent_lib.fasta', Acontent_lib)



#Final_lib different than strand with lower A content for some reason, 
#difference is in primers, lower A content file is correct and was sent to 
#sequencing

#Read in cre_test_final_lib, remove primers and add actual primers

old_primer_lib = csv_reader('cre_test_cleaned_lib.csv')

trimmed_lib = {}

for x in old_primer_lib:
    sequence = old_primer_lib[x]
    trimmed_seq = sequence[26:-26]
    
    trimmed_lib[x] = trimmed_seq
    
final_lib_update = {}
    
for x in trimmed_lib:
    sequence = trimmed_lib[x]
    final_seq_update = fwd_primer + MluI + sequence + KpnI + rev_primer
    
    final_lib_update[x] = final_seq_update
    
    csv_writer(final_lib_update, 'cre_test_final_lib_update.csv')
    



           
#oldcode

#------------------------------------------------------------------------------
#Design an optimized library that tests 5 sites with spacing that places CREs
#either at one helical turn upstream or 2 helical turns. This results in 
#either alternating 2-3 bp spacing between the 8 bp CRE sites or 9 bp spacing
#between 12 bp CRE sites with flanks. From 5' to 3' site 1 and site 5 will just
#be consensus, 2 and 3 will be 1 of the 4 affinities and 4 will be either 
#moderate or consensus. Then each arrangement will then be moved in 1 bp
#increments to the 5' end of the oligo for 10 bp total movement. This will be
#performed on a subset of 3 backgrounds. 

site_2_3 = [nosite, weak, moderate, consensus]
site_2_3_names = ['nosite', 'weak', 'moderate', 'consensus']

site_4 = [moderate, consensus]
site_4_names = ['moderate', 'consensus']

spacing_2bp = {}

for x in optim_enhancers:
    enh = optim_enhancers[x]
    for j in range(len(site_2_3)):
        for k in range(len(site_2_3)):
            for l in range(len(site_4)):
                
                site1 = consensus
                site2 = site_2_3[j]
                site3 = site_2_3[k]
                site4 = site_4[l]
                site5 = consensus
                
                seqs = [enh[:i] + 'AT' + site1 + 'GCT' + site2 + 'GT' + site3 
                        + 'GAT' + site4 + 'GT' + site5 + 'GC' + enh[i+54:]
                for i in range(len(enh)- 54 - 10 + 1, len(enh)- 54 + 1, 1)]
                
                #loc = seqs.find('xxxxxxxx')
                #while loc > -1:
                    #seqs = seqs[:loc] + seqs[loc:loc+8] + seqs[loc+8]
                    #loc = seq.find('xxxxxxxx')
                
                spacing_name = {'optim_turn_1_consensus_' + names[j] + '_' +
                                names[k] + '_' + names[l] + '_consensus_dist_' 
                                + str(11 - i) + '_' + x : seqs[i] 
                                for i in range(len(seqs))}
                
                spacing_2bp.update(spacing_name)
                subpool_libraries['minP'] = spacing_2bp
                subpool_libraries['shortminP'] = spacing_2bp
                print(len(spacing_2bp), "optimized 1 turn sequences in both subpools")
    
csv_writer(spacing_2bp, 'cre_test_optim_1.csv')









def name_gc_content_background(backgrounds):
    gc_names = []
    name_number = 1
    seq = list(backgrounds.keys())
    seq_name = {}
    
    for gc in backgrounds.values():
        gc_name = 'background_' + str(name_number) + '_gc_' + str(gc)
        name_number = name_number + 1
        gc_names.append(gc_name)
        
    for i in range(len(seq)):
        seq_name[gc_names[i]] = seq[i]
        
    return(seq_name)
    
bgs_top4_name = name_gc_content_background(bgs_top4)


#This design will be present in both subpools to additionally test the effect 
#of distance to a shorter minimal promoter

















