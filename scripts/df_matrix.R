#!/usr/bin/env Rscript

#Taking in arguments
library("optparse")
 
option_list = list(
  make_option(c("-f", "--filenames"), type="character", default=NULL, 
              help="order in which user prefers to align the matrix",metavar="character"),
  make_option(c("-c", "--colnum"), type="integer", default=4, 
              help="order in which user prefers to align the matrix",metavar="integer"),
  make_option(c("-s", "--strandcol"), action = "store_true", default = FALSE,
              help="will take 6 columns instead of 3 columns of bed file for merging",metavar="logical"),
  make_option(c("-a", "--all"),  action = "store_true", default = FALSE, 
              help="keep all entries, file may be large based on number of samples being used",metavar="logical"),
  make_option(c("-t", "--thresh"), type="numeric", default=0.2, 
              help="fraction of samples containing ribos will be filtered in the output, only works with -a flag",metavar="numeric"),
  make_option(c("-o", "--output_prefix"), type="character", default="out", 
              help="output file name [default %default]", metavar="character")
);
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
thresh=opt$thresh

#Checking for required input
writeLines("\n...Checking input...\n")
if (is.null(opt$filenames)){
  print_help(opt_parser)
  stop("Please specify annotated file.n", call.=FALSE)
}
if (is.null(opt$colnum)){
  print_help(opt_parser)
  stop("Please specify the column to be used to make matrix.n", call.=FALSE)
}


#####Main Script
#1.Calling required Packages
writeLines("\n...Calling required Package(s)...\n")

library(GenomicRanges, quiet = T) 
library(plyranges, quiet = T)

#2.Defining functions if any
writeLines("\n...Defining required functions...\n")

#3. Preprocessing input files
writeLines("\n...Processing Input files...\n")

filenames = read.table(opt$filenames, header=F, sep="\t", quote="")
filepatterns=filenames[,"V1"]


#4.Main code
writeLines("\n...Executing Main code...\n")

for (p in filepatterns) {
  writeLines(paste("Working on", p))
  seq.df=read.table(list.files(pattern=p)[1], sep='\t', header=F)
  colnames(seq.df)[1:3]<-c("chr","start","end")
  colnames(seq.df)[opt$colnum]=p

  if (opt$strandcol) { 
    df=seq.df[,c(1:6,opt$colnum)]
  } else {
    df=seq.df[,c(1:3,opt$colnum)]
  }
  
  if (p == filepatterns[1]){
    mat=df
  } else {
    if (opt$all) {
      mat=merge(mat,df, all=TRUE)
    } else {
      mat=merge(mat,df)
    }
    
  }
}

if (opt$strandcol) { 
    mat$samplexists<-apply(mat[7:ncol(mat)], 1, function(x) sum(!is.na(x)))
    if (thresh > 1) {
      mat_filtered<-mat[mat$samplexists>=thresh,]
    } else {
      mat_filtered<-mat[mat$samplexists>=thresh*(ncol(mat)-7),]
    }
  } else {
    mat$samplexists<-apply(mat[4:ncol(mat)], 1, function(x) sum(!is.na(x)))
    if (thresh > 1) {
      mat_filtered<-mat[mat$samplexists>=thresh,]
    } else {
      mat_filtered<-mat[mat$samplexists>=thresh*(ncol(mat)-4),]
    }
  }

write.table(mat_filtered, file=opt$output_prefix, sep = "\t",quote=F,row.names=F, col.names=T)






#5.Plotting

writeLines("\nJob Finished.\n")
writeLines(paste("Output file name:", opt$output_prefix))
