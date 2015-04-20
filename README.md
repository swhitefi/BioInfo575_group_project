# BioInfo575 group project: Phylogenetic methods used in Bioinfo 575 final project Spring 2015

##Core Genome
###### download parsnp and gingr
	curl -O https://github.com/marbl/harvest-tools/releases/download/v1.2/harvesttools-Linux64-v1.2.tar.gz/
###### all documentation is available here: http://harvest.readthedocs.org/en/latest/content/harvest-tools.html
######perform alignment with forced inclusion of all genomes in directory and random reference genome
###### input is a directory of genomes 
		parsnp -r ! -d ./ -p 8 -c
###### visualized Newick formatted trees in Gingr following instructions above and R as specified below.
##Pan Genome
### download LS-BSR and documentation
### Download manual
	curl -O https://github.com/jasonsahl/LS-BSR/blob/master/LS_BSR_manual.pdf?raw=true
### follow instructions in the manual to clone the repository with git to download LS-BSR programs
	git clone https://github.com/jasonsahl/LS-BSR.git
#### Follow instructions in manual to download dependencies 
### Run LS-BSR
	python ls_bsr.py -d directory_of_FASTA -u Path/to/USEARCH
###Look at pan genome stats: klebsiella example:
	python pan_genome_stats.py -b klebs_bsr_matrix_values.txt 
### output:
	Number of conserved genes = 4007
	Number of of unique genes = 468
	Number of unique genes per genome = 23.4
### Plot LSBSR matrices in R
#####ABAU
	Abau_LSBSR_mat<-read.delim("./ABAU/Abau_bsr_matrix_values.txt", row.names=1)
	Abau_LSBSR_mat<-as.matrix(Abau_LSBSR_mat)
######find variance
	library(matrixStats)
	test<-sort(rowVars(Abau_LSBSR_mat), decreasing=TRUE)
	getvar<-apply(Abau_LSBSR_mat[,-1],1,var)
	param<-.00001
	test<-Abau_LSBSR_mat[getvar > param & !is.na(getvar), ]
#######plot
	rc<-rainbow(nrow(test), start = 0, end = .5)
	cc<-rainbow(ncol(test), start = 0, end = .5)
	heatmap(test, col=cm.colors(100),scale="column",RowSideColors=rc, ColSideColors = cc, margins=c(11,11), main="Abau 		LSBSR",cexRow=.3, cexCol=.8)
#######KLEBS
	Klebs_LSBSR_mat<-read.delim("./KLEBS/klebs_bsr_matrix_values.txt", row.names=1)
	Klebs_LSBSR_mat<-as.matrix(Klebs_LSBSR_mat)
	red_Klebs_LSBSR_mat<-Klebs_LSBSR_mat[1:100,]
	plot(sort(red_Klebs_LSBSR_mat))
	rc<-rainbow(nrow(red_Klebs_LSBSR_mat), start = 0, end = .5)
	cc<-rainbow(ncol(red_Klebs_LSBSR_mat), start = 0, end = .5)
	heatmap(red_Klebs_LSBSR_mat, col=cm.colors(100),scale="column",RowSideColors=rc, ColSideColors = cc, margins=c(11,11), 		main="Klebs LSBSR",cexRow=.3, cexCol=.5)
#######select ones with highest variance
	getvar<-apply(Klebs_LSBSR_mat[,-1],1,var)
	param<-.01
	test<-Klebs_LSBSR_mat[getvar > param & !is.na(getvar), ]
	rc<-rainbow(nrow(test), start = 0, end = .5)
	cc<-rainbow(ncol(test), start = 0, end = .5)
	heatmap(test, col=cm.colors(100),scale="column",RowSideColors=rc, ColSideColors = cc, margins=c(11,11), main="Klebs 		LSBSR",cexRow=.05, cexCol=.5)

##Kmer Based SNP and kmer Based SNP Core Genome
######download kSNP3 and dependencies
	curl -O http://sourceforge.net/projects/ksnp/files/kSNP3.0_Linux_package.zip/download
#######downloaded kSNP3 to flux

#######then put all genomes from genome directory into a single fasta

#######the kchooser program requires a fasta with all genomes where each sequences is on one line

#######run the program (included) to remove the lines
	/path/kSNP3_Linux_package/kSNP3/fasta_remove_new_lines3 
	/path/Abau_all_genomes.txt > Abau_all_nolines.fasta
####### the output is a file called Abau_all_nolines.fasta
	/home2/swhitefi/Bio_info.packages/kSNP3_Linux_package/kSNP3/fasta_remove_new_lines					 	/home2/swhitefi/Bio_info.packages/Genomes/KLEBS_GENOMES/Kleb_ALL_8-20-12_anon/Klebs_all_genomes.txt 				 Klebs_all_nolines.fasta
####### the kchooser program reported that:
	the optimum value for kmer length for Klebs is 25
	the optimum value for kmer length for ABAU is 19
	flags:
	-in inputfile 	file listing paths and genome names of genome fasta files for analysis
	-outdir name 	of directory to put outfiles in
	-k 				integer (odd) length of k-mers flanking SNPS
	-vcf 			output a vcf file
	-ML 			calculate a Maximum likelihood tree
	-annotate 		to get SNP positional info. needs a file containing genome names
	(also run with -core to calculate tree based only on SNPS in core genome)
	- | tee Run1Log  to get logfile

#######first make infile with script that comes with kSNP so the infile is tab delimited
####### with path and name of the genome
#######go to directory with the genomes in it then:
	MakeKSNP3infile ./ Abau_infile A
	MakeKSNP3infile ./  klebs_infile A
#######run kSNP3:
	/path/kSNP3_Linux_package/kSNP3/kSNP3 
	-in /path/kSNP3_output/klebs_output/Abau_infile 
	-k 19 -outdir /path/Bio_info.packages/kSNP3_output/klebs_output/kSNP_klebs_output 
	-vcf -ML
######download figtree to visualize the trees!
	-curl -O http://tree.bio.ed.ac.uk/download.php?id=90
##SNP phylogeny
######download mugsy and setup dependencies
	curl -O https://sourceforge.net/projects/mugsy/files/
###### run mugsy. mugsy needs all genomes and paths called in command line. 
	mugsy --directory ./mugsy_output/ -p Abau_test ~/Bio_info.packages/Genomes/ABAU_GENOMES/Abau_ALL/AbauA_genome.fasta 		~/Bio_info.packages/Genomes/ABAU_GENOMES/Abau_ALL/AbauB_genome.fasta 								~/Bio_info.packages/Genomes/ABAU_GENOMES/Abau_ALL/AbauC_genome.fast								 ~/Bio_info.packages/Genomes/ABAU_GENOMES/Abau_ALL/ACICU_genome.fast								 ~/Bio_info.packages/Genomes/ABAU_GENOMES/Abau_ALL/Abau_HC64_021406_genome.fasta
###### if lots of genomes run this python script to generate mugsy command line calls
	#!/usr/bin/env python
	from __future__ import print_function
	import argparse
	import os
	import re
	import sys
	parser = argparse.ArgumentParser(description='''This program takes a directory of genomes
	and the name of a .PBS file and creates a PBS file to run Mugsy to perform multiple genome
	alignment on all of the .fasta genome files in the directory.  After you run this script
	you will need to specify the path to mugsy in the output file
	and you will want to change the name of your flux run in the PBS preamble''')
	parser.add_argument('-GD', help='Path to genome directory')
	parser.add_argument('-mOD', help = 'Path to mugsy output directory')
	parser.add_argument('-p', help = 'prefix to name .MAF mugsy output file')
	parser.add_argument('-o', '--outfile', help='''Name of .PBS file to make.''')
	args=parser.parse_args()
	#get genome files from the directory specified in -GD
	genome_files = os.listdir(args.GD)
	outfile = open(args.outfile, 'w')
	GD=args.GD
	#print commands to run Mugsy
	#variable for directory to output data to
	mugsyOutDir = args.mOD
	prefix = args.p
	#prefix to name .MAF file
	mugsStuff = " ".join(["mugsy --directory", mugsyOutDir, "-p", prefix])
	#print the genomes
	for i in genome_files:
		print('/'.join([GD,i]), end=' ', file = outfile)
	outfile.close()
######Mugsy outputs a .MAF file

#######now need a script to convert MAF to FASTA to build a tree with fasttree
#######use mugsy maf2fasta.pl script
#######Extract an alignment, e.g., the first LCB (locally co-linear block):
	~/Bio_info.packages/mugsy_x86-64-v1r2.3/maf2fasta.pl 1 < ./klebs_mugsy.maf > test_klebs.fasta
#######now need to remove "= sign" from end of files
	grep -v "=" test_klebs.fasta > final_klebs.fasta
#######now build a tree with FastTree
	/path/FastTree -nt < final_klebs.fasta > Klebs_snp.tree
#######Trees were visualized with Figtree and R as described above for other methods

### Compare Parsnp core and kSNP core agreement in R

	library(ape)
	#read in the trees
	#klebsiella
	parsnp<-read.tree(file="Klebs_parsnp.tree")
	kSNP_Core<-read.tree(file="tree.core.tre")
	#Abau
	Abau_parsnp<-read.tree(file="Abau_parsnp.tree")
	Abau_kSnp_Core<-read.tree(file="Abau_tree.core.tree")

	#make distance matrices
	#klebs
	parsnp_mat<-cophenetic(parsnp)
	kSNP_core_mat<-cophenetic(kSNP_Core)
	#abau
	Abau_parsnp_mat<-cophenetic(Abau_parsnp)
	Abau_kSNP_core_mat<-cophenetic(Abau_kSnp_Core)
	#sort the matrices to get them in the same order
	#klebs
	parsnp_mat<-parsnp_mat[order(colnames(parsnp_mat)),order(rownames(parsnp_mat))]
	kSNP_core_mat<-kSNP_core_mat[order(colnames(kSNP_core_mat)),order(rownames(kSNP_core_mat))]
	#Abau
	Abau_parsnp_mat<-Abau_parsnp_mat[order(colnames(Abau_parsnp_mat)),order(rownames(Abau_parsnp_mat))]Abau_kSNP_core_mat<-		Abau_kSNP_core_mat[order(colnames(Abau_kSNP_core_mat)),order(rownames(Abau_kSNP_core_mat))]
	#Klebsiella & Abau: PLOT KSNP CORE, PARSNP
	plot(kSNP_core_mat[lower.tri(kSNP_core_mat)],(parsnp_mat[lower.tri(parsnp_mat)]), type="p", col="black", 			main="Correlation Between Core Genome Methods", ylab="Parsnp Core Genomic Distance", xlab="kSNP Core Genomic 			Distance")
	points(Abau_kSNP_core_mat[lower.tri(Abau_kSNP_core_mat)],(Abau_parsnp_mat[lower.tri(Abau_parsnp_mat)]), type="p",pch=4, 	col="green")
	#green = Abay
	#black = klebsiella
	#correlation statistic between matrices
	cor.test(kSNP_core_mat,parsnp_mat)
	cor.test(Abau_kSNP_core_mat,Abau_parsnp_mat)

##Bayesian Inference
##### MrBayes accepts nexus aligned sequence files as input. Therefore use Mugsy and Maf2Fasta like we did for the SNP trees to align the sequences. And then follow the next few steps to convert the fasta files into nexus files.
###### Download bioscripts.convert 4
	curl -O https://pypi.python.org/packages/source/b/bioscripts.convert/bioscripts.convert-0.4.tar.gz
	tar zxvf bioscripts.convert.tgz
	cd bioscripts.convert-0.4
	python setup.py install 
###### To convert the fasta files into nexus files
	cd bioscripts
	cd convert 
	python convbioseq.py nexus <input file>
###### Download and install MrBayes (version 3.2.4 x 64)
	curl -O http://mrbayes.sourceforge.net/download.php
###### Go into the directory where mrbayes is then ...
	./configure
	make
	make install
###### Now run Mrbayes using the following commands
	mb 
	execute <nexus file location>
	lset nst=6 rates=invgamma
	mcmc ngen=20000 samplefreq=100 printfreq=100 diagnfreq=1000
	while standard deviation > 0.01:
		yes
		10000
	no
	sump 
	sumt
###### The resulting nexus consensus trees can be visualized using Figtree or in R using the ape package
	
## R-Code with Ape package for visualization of Acinetobacter baumannii trees
	library(ape)
###### Read in tree files
	Abau_kmer <- read.tree(file)
	Abau_snp <- read.tree(file)
	Abau_core <- read.tree(file)
	Abau_BI <- read.nexus(file)
###### Plot Kmer tree with scale bar 
	plot(Abau_kmer)
	title ('Kmer Tree')
	locator()
	add.scale.bar(x=1.2, y=1)
###### Plot SNP tree with scale bar
	plot(Abau_snp)
	title('SNP Tree')
	locator()
	add.scale.bar(x=0.001, y=1)
###### Plot Core Genome tree with scale bar
	plot (Abau_core)
	title('Core Tree')
	add.scale.bar ()
###### Plot Bayesian Inference tree with scale bar
	plot(Abau_BI)
	title('Bayesian Inference Tree')
	add.scale.bar(x=0.0005, y=1.15)
###### Useful function to find location for scale bar
	locator()
## Figtree was used to visualize the Klebsiella pneumoniae trees
##### Figtree was downloaded from http://tree.bio.ed.ac.uk/software/figtree/
###### The trees were then opened in Figtree and the branches were transformed to be proportional
###### The line width was changed to 2 and a scale bar was added.
##### Photoshop was used to combine the four trees, add titles, and make, the labels easier to read.
