# BioInfo575_group_project: Phylogenetic methods used in Bioinfo 575 final project Spring 2015

##MCMC
Komal type here:
###instruction
#### do another thing
##### code comments
	your code here something


##Core Genome

##Pan Genome
### dowload LS-BSR and documentation
### Downloaded manual
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




