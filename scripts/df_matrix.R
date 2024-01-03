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
  make_option(c("-t", "--thresh"), type="numeric", default=0.1, 
              help="fraction of samples containing ribos will be filtered in the output, only works with -a flag",metavar="numeric"),
  make_option(c("-o", "--output_prefix"), type="character", default="out", 
              help="output file name [default %default]", metavar="character")
);
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
file=opt$filenames; print(file)
col=opt$colnum; print(col)
all=opt$all; print(all)
thresh=opt$thresh; print(thresh)
strand=opt$strandcol; print(strand)
out=opt$output_prefix; print(out)

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

filenames = read.table(file, header=F, sep="\t", quote="")
filepatterns=filenames[,"V1"]


#4.Main code
writeLines("\n...Executing Main code...\n")

#Merge matrix
for (p in filepatterns) {
  writeLines(paste("Working on", p))
  seq.df=read.table(list.files(pattern=p)[1], sep='\t', header=F)
  colnames(seq.df)[1:3]<-c("chr","start","end")
  colnames(seq.df)[col]=p

  if (strand) { 
    df=seq.df[,c(1:6,col)]
  } else {
    df=seq.df[,c(1:3,col)]
  }
  
  if (p == filepatterns[1]) {
    mat=df
  } else {
    if (all) {
      mat=merge(mat,df, all=TRUE)
    } else {
      mat=merge(mat,df)
    }

  }

}

nrow(mat)

#Counts the number for samples having values and apply threshold for filtration
mat=mat[mat$chr!= '2micron',]
nrow(mat)

n=ncol(mat)
if (strand) {
  mat$samplexists<-apply(mat[7:n], 1, function(x) sum(!is.na(x)))
  ncol(mat)
  mat$average<-rowMeans(replace(mat[7:n], is.na(mat[7:n]), 0))
  ncol(mat)
  mat=mat[order(mat$average ,decreasing = TRUE),]
  mat=mat[order(mat$samplexists ,decreasing = TRUE),]
  if (thresh > 1) {
    mat_filtered<-mat[mat$samplexists>=thresh,]
  } else {
    mat_filtered<-mat[mat$samplexists>=thresh*(n-8),] #thresh=1 will give rNMPs present in all samples
    ncol(mat_filtered)
  }
  } else {
  mat$samplexists<-apply(mat[4:n], 1, function(x) sum(!is.na(x)))
  mat$average<-rowMeans(replace(mat[4:n], is.na(mat[4:n]), 0))
  mat=mat[order(mat$average ,decreasing = TRUE),]
  mat=mat[order(mat$samplexists ,decreasing = TRUE),]
  if (thresh > 1) {
    mat_filtered<-mat[mat$samplexists>=thresh,]
  } else {
    mat_filtered<-mat[mat$samplexists>=thresh*(ncol(mat)-5),]
  }
}



#5.Output
#writeLines("\n...Distribution of average counts or EF...\n")

#table(mat$average)
writeLines(file)
writeLines("\n...Composition for mentioned threshold...\n")
table(mat_filtered$V4)

writeLines("\nJob Finished.\n")
writeLines(paste("Output file name:", out))
write.table(mat, file=paste(out,".all", sep=""), sep = "\t",quote=F,row.names=F, col.names=T)
write.table(mat_filtered, file=out, sep = "\t",quote=F,row.names=F, col.names=T)
