#!/usr/bin/env Rscript

#Taking in arguments
library("optparse")
 
option_list = list(
  make_option(c("-r", "--ribo"), type="character", default=NULL, 
              help="ribo bed file with col: chr start stop . . strand",metavar="character"),
  make_option(c("-c", "--cov"), type="character", default=NULL, 
              help="DNAseq coverage file with 3 col: chr location count",metavar="character"),
  make_option(c("-g", "--gen"), type="character", default=NULL, help="genome size(.fai) file", metavar="character"),
  make_option(c("-o", "--output_prefix"), type="character", default="out", 
              help="output file name [default %default]", metavar="character")
);
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#Checking for required input
writeLines("\n...Checking input...\n")
if (is.null(opt$ribo)){
  print_help(opt_parser)
  stop("Please specify ribo bed file.n", call.=FALSE)
}
if (is.null(opt$cov)){
  print_help(opt_parser)
  stop("Please specify coverage file.n", call.=FALSE)
}

if (is.null(opt$gen)){
  print_help(opt_parser)
  stop("Please specify genome(.fai) file.n", call.=FALSE)
}


ribo=opt$ribo; writeLines(ribo)
cov=opt$cov; writeLines(cov) #cov='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/AGS/DNAseq/aligned/Y1.cov'
gen=opt$gen; writeLines(gen) #gen='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/sacCer2/sacCer2-nucl.fa.fai'
out=opt$output_prefix; writeLines(out)

#####Main Script
#1.Calling required Packages
writeLines("\n...Calling required Package(s)...\n")

suppressMessages(library(GenomicRanges, quiet = T))
suppressMessages(library(plyranges, quiet = T))
suppressMessages(library(data.table, quietly = T))
suppressMessages(library(tools,quietly = T)) #for getting file basename

pre=file_path_sans_ext(basename(ribo))


#2.Defining functions if any
#writeLines("\n...Defining required functions...\n")

#3. Preprocessing input files
#writeLines("\n...Processing Input files...\n")
genome =  as.list(read.table(gen, header=F, sep="\t", quote="")[1])$V1
ribo.df = read.table(ribo, header=F, sep="\t", quote=""); ribo.df=ribo.df[,c(1,3,6)]
ribo.df=data.frame(setDT(ribo.df)[,list(count=.N),names(ribo.df)]); colnames(ribo.df)<-c("chr","loc","strand","counts") ; ribo.df = ribo.df[ribo.df$chr %in% genome,] #Getting counts data for ribo bed file
cov.df = read.table(cov, header=F, sep="\t", quote=""); colnames(cov.df)<-c("chr","loc","cov") ; cov.df = cov.df[cov.df$chr %in% genome,]; cov.df$cov<-cov.df$cov/mean(cov.df$cov) #Getting relative coverage for each genome base

#4.Main code
#writeLines("\n...Executing Main code...\n")
final.df=merge(ribo.df,cov.df,by=c("chr","loc")); if (nrow(ribo.df)==nrow(final.df)) {writeLines("Merge of ribos and coverage sucessfull")}
min_cov=min(final.df$cov[final.df$cov>0]); final.df$cov[final.df$cov==0]<-min_cov #coverage correction
final.df$norm_count=final.df$counts/final.df$cov
final.df$round_norm_counts<-round(final.df$norm_count) #Getting normalized counts
avg_ribo=sum(final.df$norm_count)/sum(read.table(gen, header=F, sep="\t", quote="")[2]) #Getting normalized counts
final.df$EF<-final.df$norm_count/avg_ribo
final.df$start<-final.df$loc - 1
final.df=final.df[,c("chr", "start", "loc", "counts", "round_norm_counts", "strand", "norm_count", "EF")]
writeLines(paste("Percent of ribo locations changed",nrow(final.df[final.df$counts != final.df$round_norm_counts,])/nrow(final.df)*100,"%")) 
writeLines(paste("Percent of ribo locations tend to zero",nrow(final.df[final.df$round_norm_counts==0,])/nrow(final.df)*100,"%")) 

#writeLines("\nJob Finished.\n")
#writeLines(paste("Output folder name:", out))
#writeLines(paste("Output columns:", c("chr", "start", "loc", "counts", "round_norm_counts", "strand", "norm_count", "EF")))
write.table(final.df, file=paste(out,"/",pre,".norm.counts", sep=""), sep = "\t",quote=F,row.names=F, col.names=F)
write.table(final.df[final.df$round_norm_counts>0,], file=paste(out,"/",pre,".bed", sep=""), sep = "\t",quote=F,row.names=F, col.names=F)
