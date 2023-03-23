#!/usr/bin/env Rscript

library("optparse")

option_list = list(
  make_option(c("-f", "--bed_file"), type="character", default=NULL, 
              help="ranges to be annotated in bed format with headers for first three columns as follows: chr\tstart\tend \n ", metavar="character"),
  make_option(c("-o", "--output_folder"), type="character", default="out", 
              help="output file name [default %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

writeLines("\n...Checking input...\n")


writeLines("Loading packages")
library(data.table, quiet = T)
library(dplyr, quiet = T)
library(stringr, quiet = T)
#library(org.Hs.eg.db, quiet = T)
library(annotatr, quiet = T)
library(stringr, quiet = T)
library(GenomicRanges, quiet = T) 
library(biomaRt, quiet = T)
library(Biostrings, quiet = T)
library(plyranges, quiet = T)
suppressMessages(library(tools,quietly = T)) #for getting file basename
chrord<-c(paste("chr", seq(1:22), sep=""), "chrX", "chrY")

writeLines("Creating annotations")
##### Cpg annotations
cpg_elements=sort(annotatr::build_annotations(genome = "hg38", annotations = "hg38_cpgs"))
cpg_elements.anno=keepSeqlevels(cpg_elements,chrord, pruning.mode="coarse")
cpg_elements.anno=cpg_elements.anno[,c("type")]#<<<<<<<<<<<<<<<<<<< used for final annotation

##### Genic element annotations
#Ens104.trans.db<-read.table('/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/hu-anno/Ensemble/ENS_Homo_sapiens.GRCh38.103.transcripts.bed', sep="\t", header=F)
#colnames(Ens104.trans.db)<-c("chr", "start", "end", "strand","tx_id", "gene_id", "tx_name")
#Ens104.genes.db<-read.table('/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/hu-anno/Ensemble/ENS_Homo_sapiens.GRCh38.104.genes.bed', sep="\t", header=F)
#colnames(Ens104.genes.db)<-c("chr", "start", "end", "strand","gene_id", "Gene_name", "type", "subtype", "description")

#Genenameanno<-merge(Ens104.trans.db[,c("tx_id", "gene_id")],
      #Ens104.genes.db[,c("gene_id", "Gene_name", "type", "subtype", "description")],
      #by="gene_id")

#gene_elements=sort(annotatr::build_annotations(genome = 'hg38' , annotations = 'hg38_basicgenes'))
#mcols(gene_elements)$str<-strand(gene_elements)
#mcols(gene_elements)$tx_id<-sapply(strsplit(mcols(gene_elements)$tx_id, "\\."), "[[", 1)
#gene_elements=gene_elements[,c("type","tx_id", "str")]
#gene_elements.df<-data.frame(gene_elements)
#gene_elements.df=merge(gene_elements.df,Genenameanno, by="tx_id" )
#gene_elements.anno=sort(makeGRangesFromDataFrame(gene_elements.df[,-1],keep.extra.columns=TRUE))#<<<<<<<<<<<<<<<<<<< used for final annotation

gene_elements.db<-read.table('/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/anno/standardanno/combined_genic_polyA_anno.bed', sep="\t", header=F, quote="")
#gene_elements.db<-gene_elements.db[,c(1:7,9)] #Don't need more right now
gene_elements.db<-gene_elements.db[,c(1:3,6,7)]
#colnames(gene_elements.db)<-c("chr", "start", "end", "length", "ENST ID","strand","generegion", "genename")
colnames(gene_elements.db)<-c("chr", "start", "end","strand","generegion")
gene_elements.anno<-makeGRangesFromDataFrame(gene_elements.db,keep.extra.columns=TRUE, starts.in.df.are.0based = TRUE)#<<<<<<<<<<<<<<<<<<< used for final annotation

##### CCRE regions
ccre.db<-read.table('/storage/coda1/p-fstorici3/0/dkundnani3/rich_project_bio-storici/Hu_analysis/anno/standardanno/combined_ccre_annno.bed', sep="\t", header=F)
ccre.db<-ccre.db[,c(1:3,6,11)] 
colnames(ccre.db)<-c("chr", "start", "end", "strand","CCRE")
ccre.anno<-makeGRangesFromDataFrame(ccre.db,keep.extra.columns=TRUE, starts.in.df.are.0based = TRUE)


##### Repeat Regions
rmsk.db<-read.table('/storage/coda1/p-fstorici3/0/dkundnani3/rich_project_bio-storici/hu-anno/UCSC/repeat_mask.bed', sep="\t", header=F)
#rmsk.db<-rmsk.db[,c(6:8,10:13)] 
rmsk.db<-rmsk.db[,c(6:8,10,12)]
#Don't need more right now
#colnames(rmsk.db)<-c("chr", "start", "end", "strand","repName", "repClass", "repFamily")
colnames(rmsk.db)<-c("chr", "start", "end", "strand","repClass")
rmsk.anno<-makeGRangesFromDataFrame(rmsk.db,keep.extra.columns=TRUE, starts.in.df.are.0based = TRUE)#<<<<<<<<<<<<<<<<<<< used for final annotation

##### ERIZ regions
#eriz1.db<-read.table('/storage/coda1/p-fstorici3/0/dkundnani3/rich_project_bio-storici/hu-anno/ARIZ/hglft_hg38_ERIZ_K652.bed', sep="\t", header=T)
#eriz1.gr<-makeGRangesFromDataFrame(eriz1.db,keep.extra.columns=TRUE, starts.in.df.are.0based = TRUE)#<<<<<<<<<<<<<<<<<<< used for final annotation
#eriz2.db<-read.table('/storage/coda1/p-fstorici3/0/dkundnani3/rich_project_bio-storici/hu-anno/ARIZ/hglft_hg38_ERIZ_GM12878.bed', sep="\t", header=T)
#eriz2.gr<-makeGRangesFromDataFrame(eriz2.db,keep.extra.columns=TRUE, starts.in.df.are.0based = TRUE)#<<<<<<<<<<<<<<<<<<< used for final annotation

##### ORM regions
#orm.db<-read.table('/storage/coda1/p-fstorici3/0/dkundnani3/rich_project_bio-storici/hu-anno/ORM/IZ_ORM_hg38.bed', sep="\t", header=T)
#orm.gr<-makeGRangesFromDataFrame(orm.db,keep.extra.columns=TRUE, starts.in.df.are.0based = TRUE)#<<<<<<<<<<<<<<<<<<< used for final annotation

##### Blacklist regions
blacklist.anno<-makeGRangesFromDataFrame(read.table('/storage/coda1/p-fstorici3/0/dkundnani3/rich_project_bio-storici/hu-anno/ENCODE/GRCh38_unified_blacklist.bed', sep="\t", header=T),keep.extra.columns=TRUE, starts.in.df.are.0based = TRUE) 
#<<<<<<<<<<<<<<<<<<< used for final annotation

##### Bringing it together

seq.df<-read.table(opt$bed_file,sep="\t",header=F)
colnames(seq.df)[1:3]<-c("chr","start","end")
colnames(seq.df)[6]<-"strand"
pre=file_path_sans_ext(basename(opt$bed_file))
writeLines("Annotating")
seq.gr<-sort(makeGRangesFromDataFrame(seq.df,keep.extra.columns=TRUE, starts.in.df.are.0based=TRUE))
#annoted<-unique(join_overlap_left(join_overlap_left(join_overlap_left(join_overlap_left(join_overlap_left(join_overlap_left(join_overlap_left(seq.gr, cpg_elements.anno), gene_elements.anno), rmsk.anno), eriz1.gr), eriz2.gr),orm.gr), blacklist.gr))
annoted<-join_overlap_left(join_overlap_left(join_overlap_left(join_overlap_left(join_overlap_left(seq.gr, cpg_elements.anno), gene_elements.anno), ccre.anno),rmsk.anno), blacklist.anno)
start(annoted)=start(annoted)-1
annoted.df<-distinct(data.frame(annoted))
write.table(annoted.df, paste(opt$output_folder,"/",pre,"_annotated.bed", sep=""), sep = "\t",quote=F,row.names=F, col.names=F)
writeLines(paste("Annotated file saved in", paste(opt$output_folder,"/",pre,"_annotated.bed", sep="")))


