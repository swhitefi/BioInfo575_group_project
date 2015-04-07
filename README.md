# BioInfo575_group_project
Code for phylogenetic methods used in Bioinfo 575 final project Spring 2015

Phylogenetic tree building methods

MCMC

Core Genome

Pan Genome
LS-BSR 
## Downloaded manual

	curl -O https://github.com/jasonsahl/LS-BSR/blob/master/LS_BSR_manual.pdf?raw=true

# following instructions in the manual to clone the repository with git
	
	$git clone https://github.com/jasonsahl/LS-BSR.git

# follow instructions to download dependencies
#run LSBSR
python ls_bsr.py -d directory_of_FASTA -u Path/to/USEARCH

#look at pan genome stats
#klebsiella example:
$ python pan_genome_stats.py -b klebs_bsr_matrix_values.txt 
# output:
# of conserved genes = 4007
# of unique genes = 468
# of unique genes per genome = 23.4


Kmer Based SNP and kmer Based SNP Core Genome

SNP 

