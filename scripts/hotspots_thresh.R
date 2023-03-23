#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", args[grep(file.arg.name, args)])
script.basename <- dirname(script.name)
other.name <- file.path(script.basename, "hotspot_analysis.R")
source(other.name)

library("optparse")

option_list = list(
  make_option(c("-b", "--bed_file"), type="character", default=NULL, 
              help="bed file in the format: chr\tstart\tstop\tid\tscore\tstrand \n ", metavar="filetype"),
  make_option(c("-g", "--genome"), type="character", default='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38-chromosomes.fa.fai',help="genome sizes file[default %default]", metavar="character"),
  make_option(c("-r", "--ref"), type="character", default='BSgenome.Hsapiens.UCSC.hg38', 
              help="BS genome name to be used for reference", metavar="integer"),
  make_option(c("-t", "--thresh"), type="character", default='0.001', 
              help="hotspot probability threhold,\n Recommended: 0.001 for yeast and 0.01 for human [default %default]", metavar="integer"),
  make_option(c("-o", "--output_folder"), type="character", default="out", 
              help="output file name [default %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

writeLines("\n...Checking input...\n")
if (is.null(opt$bed_file)){
  print_help(opt_parser)
  stop("Please specify input bed file.n", call.=FALSE)
} else {
  file=opt$bed_file
  writeLines(paste("Found the file", opt$bed_file))
}

file=opt$bed_file

writeLines(paste("Probability threshold value: ", opt$thresh))

if (file.exists(opt$output_folder)) {
  writeLines("Output directory exists!")
} else {
  writeLines(paste("Creating output directory", opt$output_folder))
  dir.create(opt$output_folder)
}

writeLines("\n...Calling Dependencies...\n")
suppressMessages(library(data.table, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(stringr, quietly = T))
suppressMessages(library(stringr,quietly = T))
suppressMessages(library(GenomicRanges, quietly = T) )
suppressMessages(library(plyranges, quietly = T))
suppressMessages(library(Biostrings, quietly = T))
suppressMessages(library(tools,quietly = T)) #for getting file basename

#chrord<-c(paste("chr", seq(1:22), sep=""), "chrX", "chrY")


##### Get flank sequence Function

pre=file_path_sans_ext(basename(file))


#files=list.files('.', pattern='*.bed')
writeLines("\n...Starting...\n")

writeLines("\n...Reading bed file...\n")
bed<-read.table(file,sep="\t",header=F)    

colnames(bed)<-c("chr", "start", "end", "id", "score", "strand")


writeLines("\n...Getting Flanking Sequence...\n")
#genome=BSgenome.Hsapiens.UCSC.hg38
#genome=BSgenome.Scerevisiae.UCSC.sacCer2
ref=opt$ref

seq.df=get_flank_seq(bed,3,3, ref)
seq.df=seq.df[,c(1,2,3,7,9,6,8)]


write.table(seq.df,paste(opt$output_folder,"/",pre,"_allseq.bed", sep=""), sep="\t", row.names=F, col.names=F, quote=F)

writeLines("\n...Counting...\n")
count=data.frame(setDT(seq.df)[,list(count=.N),names(seq.df)])
#counts=filter(count, chr != "chrM")
count=count[order(count$count, decreasing=TRUE),]
    
writeLines("\n...Getting distribution probabilities...\n")
N=sum(count$count)
U=length(count$count)
lamda=N/U
count$probability=ppois(count$count-1, lambda=lamda,lower=FALSE ) ### P (X > k-1) = P(X >= k)


writeLines("\n...Getting EF...\n")
genome = read.table(opt$genome, header=F, sep="\t")
genomesize=sum(genome$V2)
count$EF=(count$count/sum(count$count))*genomesize
tail(count)
#count$EF_pval=ppois(count$count-1, lambda=1,lower=FALSE )

write.table(count,paste(opt$output_folder,"/",pre,"_poissprob.bed", sep=""), sep="\t", row.names=F, col.names=F, quote=F)

writeLines("\n...Getting hotspots...\n")
##### Get hotspots using probability threshold
t=as.double(opt$thresh)
count_thresh=tail(filter(count, probability<t),1)$count-1
percent_thresh=count[round(nrow(count))*t,]$count-1

write.table(filter(count, count>percent_thresh),paste(opt$output_folder,"/",pre,"_percentthresh_",opt$thresh,".bed", sep=""), sep="\t", row.names=F, col.names=F, quote=F) 
ggsave(paste(opt$output_folder,"/out/percent/",pre,"_percentthresh_",opt$thresh,".jpeg", sep=""), plot = plot_meme(filter(count, count>percent_thresh),7))

write.table(filter(count, count>count_thresh),paste(opt$output_folder,"/",pre,"_ppthresh_",opt$thresh,".bed", sep=""), sep="\t", row.names=F, col.names=F, quote=F) 
ggsave(paste(opt$output_folder,"/out/ppthresh/",pre,"_ppthresh_",opt$thresh,".jpeg", sep=""), plot = plot_meme(filter(count, count>percent_thresh),7))

writeLines(paste(pre ,"done. Output saved in", opt$output_folder, sep=""))