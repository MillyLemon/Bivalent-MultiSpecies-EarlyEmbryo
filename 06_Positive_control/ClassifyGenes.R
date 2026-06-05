###################################################
## USER PARAMETERS (SET THESE) ####################
###################################################

mapped_file_dir="./"
mapped_file_names=c("./GSM691591_ES.H3K27.standard.bam","./GSM691593_ES.WCE.standard.bam")
experiment_names=c("GSM691591_ES.H3K27","GSM691593_ES.WCE")
experiment_MACS_file="./GSM691591_ES.H3K27.standard-W200-G600-FDR0.01-island.bed"
dots=TRUE
max_sim_reads=1000000
txdb_file="./TxDb_mm9_gencode_vM1.sqlite"
pval_cut=.001
fc_cut=1
promoter=3000
chr_fix=function(u){gsub("([^T]*)T*\\..*","chr\\1",u)}



###################################################
## DO NOT EDIT BEYOND THIS LINE ###################
###################################################

###########
#PREAMBLE #
###########

source("functions.R")
library(Rsamtools)
library(GenomicFeatures)
library(edgeR)
library(GenomicAlignments)
library(GenomicRanges)

###########
#LOAD   ###
###########

setwd(mapped_file_dir)
#We need to remember the library sizes
library_sizes=c()
for(i in 1:length(mapped_file_names)){
	file=mapped_file_names[i]
	label=experiment_names[i]
	x=readGAlignments(file)
	library_sizes=c(library_sizes,length(x))
	starts=start(x)
	chrs=rname(x)
	strands=strand(x)
	marks=paste(chrs,starts,strands,sep="_")
	tmp=table(marks)
	w=rep(match(names(tmp),marks),pmin(max_sim_reads,tmp))
	#Cleanup
	rm(marks,tmp)
	message("Filtered PCR artifacts...")
	if(dots){
		#If we're just taking a dot as representative, make sure it's at the 5' end...
		b=start(x[w])
		s=as.character(strand(x[w]))
		b[s=="-"]=end(x[w])[s=="-"]
		reads=GRanges(seqnames=rname(x[w]),ranges=IRanges(start=b,end=b),strand=strand(x[w]))
	}else{
		reads=grg(x[w])
	}
	message("Processed file ",file)
	save(reads,file=paste(label,"_reads.RData",sep=''))
}

#############
#ANNOTATION #
#############

#txdb=loadFeatures(txdb_file)
#txdb<-txdb_file
txdb <- loadDb(txdb_file)  #2025.12.5
#Calculate every transcript
tmp=unfactor(as.list(txdb))
annot=GRanges(seqnames=tmp$transcripts$tx_chrom,ranges=IRanges(start=tmp$transcripts$tx_start,end=tmp$transcripts$tx_end),strand=tmp$transcripts$tx_strand,tx_name=tmp$transcripts$tx_name,tx_id=tmp$transcripts$tx_id,gene_id=tmp$genes$gene_id[match(tmp$transcripts$tx_id,tmp$genes$tx_id)])
g2tx=values(annot)
#Now just keep the longest transcript, that's what we want
tmp=width(annot)
names(tmp)=values(annot)$tx_name
tmp=split(tmp,values(annot)$gene_id)
tmp=lapply(tmp,function(u){u[which.max(u)]})
tmp=sapply(tmp,names)
annot=annot[match(tmp,values(annot)$tx_name)]
#Save it for later use
save(annot,file="annotation.RData")


#############
#MAKE TOC   #
#############

load("annotation.RData")
#Now we need to set the extra upstream region
#The resize command anchors at the wrong end (ie, extends at the end of annotations), so we first flip the strands
strand(annot)=ifelse(as.character(strand(annot))=="+","-","+")
#Extend to include the promoter
annot=resize(annot,width(annot)+promoter)
#Flip the strands back to the proper orientation
strand(annot)=ifelse(as.character(strand(annot))=="+","-","+")
#Now make the table
region_counts=matrix(0,nrow=length(annot),ncol=length(experiment_names))
colnames(region_counts)=c(experiment_names)
rownames(region_counts)=values(annot)$gene_id
strand(annot)="*"
#Now we have to make a table of counts for this region we've defined
for(i in 1:ncol(region_counts)){
	label=colnames(region_counts)[i]
	#Need to load the relevant object
	load(paste(label,"_reads.RData",sep=''))
	strand(reads)="*"
	reads=GRanges(seqnames=chr_fix(as.character(seqnames(reads))),ranges=ranges(reads),strand=as.character(strand(reads)))
	region_counts[,i]=countOverlaps(annot,reads)
	message("Calculated overlapping reads for condition ",colnames(region_counts)[i])
}
save(region_counts,file="Poisson_test_toc.RData")

#############
#CALL GENES #
#############

#Initialize a giant table with all the genes and their calls
load("annotation.RData")
gene_classification=data.frame(gene_id=values(annot)$gene_id,Poisson=FALSE,MACS=FALSE,Classification=NA)
#For each gene, call it as bound or not by the Poisson test (and by MACS if given)
#Call using the Poisson method first
load("Poisson_test_toc.RData")
DGE=DGEList(region_counts,group=colnames(region_counts),lib.size=library_sizes)
d=estimateCommonDisp(DGE)
pvals=exactTest(d,colnames(region_counts)[2:1],dispersion=0.2) # https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/exactTest
pvals=pvals$table
pvals$FDR=p.adjust(pvals$PValue,method="BH")
poisson.bound=rownames(pvals)[pvals$FDR<pval_cut & pvals$logFC>fc_cut]
gene_classification$Poisson[gene_classification$gene_id%in%poisson.bound]=TRUE

#Now if we have MACS data, call using that too
if(length(experiment_MACS_file)!=0){
	MACS=read.table(experiment_MACS_file,sep='\t',skip=0 , header = T)
	MACS[,1]=chr_fix(MACS[,1])
	load("annotation.RData")
	#Now we need to set the extra upstream region
	#The resize command anchors at the wrong end (ie, extends at the end of annotations), so we first flip the strands
	strand(annot)=ifelse(as.character(strand(annot))=="+","-","+")
	#Extend to include the promoter
	annot=resize(annot,width(annot)+promoter)
	#Flip the strands back to the proper orientation
	strand(annot)=ifelse(as.character(strand(annot))=="+","-","+")
	#Make the MACS data into a GRANGES object
	MM=GRanges(seqnames=MACS[,1],ranges=IRanges(start=MACS[,2],end=MACS[,3]),strand="*")
	#Now find any overlaps
	set=which(countOverlaps(annot,MM)!=0)
	macs.bound=values(annot)$gene_id[set]
	gene_classification$MACS[gene_classification$gene_id%in%macs.bound]=TRUE
}

#############
#CLASIFY    #
#############

#Finally, we run the classifaction on all genes that are called.  As the classify_gene function can fail, we need to be careful
#Make a coverage object from the data...
load(paste(experiment_names[1],"_reads.RData",sep=''))
reads=GRanges(seqnames=chr_fix(as.character(seqnames(reads))),ranges=ranges(reads),strand=as.character(strand(reads)))
cover=coverage(reads)
set=which(gene_classification$MACS | gene_classification$Poisson)
stop=length(set)
stats=as.list(rep(NA,stop))
names(stats)=values(annot)$gene_id[set]
passed=rep("Unknown",stop)
genes=values(annot)$gene_id
while(length(set)!=0){
	for(j in 1:length(set)){
		i=set[j]
		tmp=tryCatch(classify_gene(genes[i],annot,cover),error=function(u){message(u);cat("\n");"FAILED"})
		stats[[match(genes[i],names(stats))]]=tmp
		if(length(tmp)==1){
			passed[j]="FAILED"
		}
		pp(length(set),j)
	}
	set=set[which(passed=="FAILED")]
	passed=rep("Unknown",length(set))
}
calls=sapply(stats,function(u){if(length(u)==1){"Unclassified"}else{u$classification}})
gene_classification$Classification=calls[match(gene_classification$gene_id,names(calls))]
#Save the classification matrix
save(stats,gene_classification,file="Gene_classification.RData")
write.table(gene_classification,row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t',file="Gene_classification.txt")


