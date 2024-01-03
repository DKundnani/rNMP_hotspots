#!/usr/bin/env Rscript


#args = commandArgs(trailingOnly=FALSE)
#file.arg.name <- "--file="
#script.name <- sub(file.arg.name, "", args[grep(file.arg.name, args)])
#script.basename <- dirname(script.name)
#other.name <- file.path(script.basename, "hotspot_analysis.R")
#source(other.name)

suppressMessages(library(optparse, quietly = T))
suppressMessages(library(ggplot2, quietly = T))
suppressMessages(library(ggseqlogo, quietly = T))
suppressMessages(library(tools,quietly = T)) #for getting file basename
suppressMessages(library(stringr,quietly = T))

option_list = list(
  make_option(c("-f", "--file"), type="character", default='FS201_nucl_percentthresh_0.01.bed', 
              help="file with flanking sequences for drawing meme \n ", metavar="character"),
  make_option(c("-c", "--column"), type="integer", default=5, 
              help="column containing the sequences [default= %default] ", metavar="integer"),
  make_option(c("-o", "--output_folder"), type="character", default="out", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

writeLines("\n...Checking input...\n")
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Please specify input file.n", call.=FALSE)
} else {
  file=opt$file
  writeLines(paste("Found the file", opt$file))
}

writeLines(paste("Column used for meme:", opt$column))

if (file.exists(opt$output_folder)) {
  writeLines("Output directory exists!")
} else {
  writeLines(paste("Creating output directory", opt$output_folder))
  dir.create(opt$output_folder)
}
file <- opt$file
col<-opt$column
pre=file_path_sans_ext(basename(file))
bed <- read.table(file,sep="\t",header=F)
colnames(bed)[1:6]<-c("chr","start","end",".",".","strand")

#bed=bed[str_length(bed[,c(opt$column)])==7,] #Make sure the length of the flank is 7 and remove those who are at the very corner of the genome. 
len=str_length(bed[,c(col)])
bed=bed[len>0,]
bed=bed[len==max(len),] #optional

if (nrow(bed)>5*10^5) {
    set.seed(12345)  
    bed = bed[sample(1:nrow(bed), 5*10^5), ]
}

#seq.df=get_flank_seq(bed,3,3)


writeLines("\n...Plotting...\n")
colnames(bed)[opt$column]<-"plotcol"

jpeg(paste(opt$output_folder,"/",pre,".jpeg", sep=""), width =1, height = 1, unit ='in', res=2400)

colors = make_col_scheme(chars = c('A', 'C', 'G', 'T'), cols = c('red2', 'blue4', 'darkorange2', 'green4'))
g=ggseqlogo(as.character(bed[,"plotcol"]), method='bits',seq_type='dna', col_scheme=colors ) + theme_classic() + xlab("") + ylab("") +
	  scale_y_continuous(limits = c(0, 2), expand = c(0.02, 0)) + scale_x_continuous(breaks = seq(1,7,1), labels = seq(-3,3,1), expand = c(0.02, 0)) +
	  theme(
          axis.title=element_blank(), axis.line = element_line(color = "black", size = 0.2), axis.ticks = element_line(color = "black", size = 0.2), axis.ticks.length = unit(0.05, "cm"), axis.text = element_text(color = "black", size = 7), plot.margin = unit(c(1,0, 0, 0), "mm")
	  )
ggsave(paste(opt$output_folder,"/",pre,".jpeg", sep=""), plot = g, scale=1, dpi=2400 )



#p=ggseqlogo(seq.df[,c("seq")], method='bits', seq_type='rna') + theme_classic() + xlab("") + ylab("") +
#	scale_x_continuous(breaks = seq(1,7,1), labels = seq(-3,3,1)) +
#	theme(
#        axis.line = element_line(color = "black", size = 1), axis.ticks = element_line(color = "black", size = 1), axis.ticks.length = unit(0.5, "cm"), axis.text = element_text(color = "black", size = 50), plot.margin = unit(c(0.5,0, 0, 0), "cm")
#	)
#ggsave(paste(opt$output_folder,"/adj_",pre,".jpeg", sep=""), plot = p)

writeLines(paste(pre,"done!"))

