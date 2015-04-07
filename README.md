# BioInfo575_group_project: Phylogenetic methods used in Bioinfo 575 final project Spring 2015

##MCMC

##Core Genome

##Pan Genome

##LS-BSR 
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
 R Abau_LSBSR_mat<-read.delim("./ABAU/Abau_bsr_matrix_values.txt", row.names=1)
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
heatmap(test, col=cm.colors(100),scale="column",RowSideColors=rc, ColSideColors = cc, margins=c(11,11), main="Abau LSBSR",cexRow=.3, cexCol=.8)
#######KLEBS
Klebs_LSBSR_mat<-read.delim("./KLEBS/klebs_bsr_matrix_values.txt", row.names=1)
Klebs_LSBSR_mat<-as.matrix(Klebs_LSBSR_mat)
red_Klebs_LSBSR_mat<-Klebs_LSBSR_mat[1:100,]
plot(sort(red_Klebs_LSBSR_mat))
rc<-rainbow(nrow(red_Klebs_LSBSR_mat), start = 0, end = .5)
cc<-rainbow(ncol(red_Klebs_LSBSR_mat), start = 0, end = .5)
heatmap(red_Klebs_LSBSR_mat, col=cm.colors(100),scale="column",RowSideColors=rc, ColSideColors = cc, margins=c(11,11), main="Klebs LSBSR",cexRow=.3, cexCol=.5)
#######select ones with highest variance
getvar<-apply(Klebs_LSBSR_mat[,-1],1,var)
param<-.01
test<-Klebs_LSBSR_mat[getvar > param & !is.na(getvar), ]
rc<-rainbow(nrow(test), start = 0, end = .5)
cc<-rainbow(ncol(test), start = 0, end = .5)
heatmap(test, col=cm.colors(100),scale="column",RowSideColors=rc, ColSideColors = cc, margins=c(11,11), main="Klebs LSBSR",cexRow=.05, cexCol=.5)
#######VRE
VRE_LSBSR_mat<-read.delim("./VRE/bsr_matrix_values.txt", row.names=1)
VRE_LSBSR_mat<-as.matrix(VRE_LSBSR_mat)
red_VRE_LSBSR_mat<-VRE_LSBSR_mat[1:100,]
plot(sort(red_VRE_LSBSR_mat))
rc<-rainbow(nrow(red_VRE_LSBSR_mat), start = 0, end = .5)
cc<-rainbow(ncol(red_VRE_LSBSR_mat), start = 0, end = .5)
heatmap(red_VRE_LSBSR_mat, col=cm.colors(100),RowSideColors=rc,scale="column", ColSideColors = cc, margins=c(11,11), main="VRE LSBSR",cexRow=.07, cexCol=.05)
#######select ones with highest variance
getvar<-apply(VRE_LSBSR_mat[,-1],1,var)
param<-.195
test<-VRE_LSBSR_mat[getvar > param & !is.na(getvar), ]
rc<-rainbow(nrow(test), start = 0, end = .5)
cc<-rainbow(ncol(test), start = 0, end = .5)
heatmap(test, col=cm.colors(100),RowSideColors=rc,scale="column", ColSideColors = cc, margins=c(11,11), main="VRE LSBSR",cexRow=.07, cexCol=.05)


##Kmer Based SNP and kmer Based SNP Core Genome

##SNP 

