#!/usr/bin/env Rscript

#Usage:source(paste(getwd(),"/hotspot_analysis.R", sep=""))
suppressMessages(library(GenomicRanges, quietly = T) )
suppressMessages(library(BSgenome, quietly = T))
#suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38, quietly = T))
#suppressMessages(library(BSgenome.Scerevisiae.UCSC.sacCer2, quietly = T))
suppressMessages(library(ggplot2, quietly = T))
suppressMessages(library(ggseqlogo, quietly = T))
suppressMessages(library(tools,quietly = T)) #for getting file basename
suppressMessages(library(stringr,quietly = T))


get_flank_seq<-function(bedfn,l,r,genome){
  bsgenome<-getBSgenome(genome)
  regions<-makeGRangesFromDataFrame(bedfn,keep.extra.columns=TRUE, seqinfo=seqinfo(bsgenome), starts.in.df.are.0based=TRUE)
  flank_regions=trim(resize(resize(regions,fix='start',l+1), fix='end', l+r+1))
  ribo=as.data.frame(BSgenome::getSeq(bsgenome, regions))
  colnames(ribo)='ribo'
  ribo=gsub("T","U", ribo[,1])
  seq=as.data.frame(BSgenome::getSeq(bsgenome, flank_regions))
  colnames(seq)='seq'
  seqU=paste(substr(seq[,1],1,l), ribo,substr(seq[,1],l+2,l+r+1), sep="")
  #ribo.rna=RNAStringSet()
  #seq.rna=list()
  #for (i in seq(1,length(seq))) {
    #if (as.character(ribo[[i]])=="T") {
      #ribo.rna[[i]]=RNAString("U")
    #} else {
      #ribo.rna[[i]]=RNAString(as.character(ribo[[i]]))
    #}
    #seq.rna[i]=paste(as.character(Views(seq[[i]], start=1:3)[[1]]),'r',as.character(ribo.rna[[i]]), as.character(Views(seq[[i]], end=7,width=3)[[1]]), sep="")
  #}
  #return(cbind(bedfn,as.data.frame(unlist(seq.rna), col.names = "Sequence"),as.data.frame(ribo.rna)))
  return(cbind(bedfn,ribo,seq, as.data.frame(seqU)))
}

plot_meme <- function(df,coln) {
  df=df[str_length(df[,c(coln)])==7,]
  if (nrow(df)>5*10^5) {
    set.seed(12345)  
    df = df[sample(1:nrow(df), 5*10^5), ]
  }
  writeLines("\n...Plotting...\n")
  colnames(df)[coln]<-"plotcol"
  #jpeg(paste(opt$output_folder,"/",pre,".jpeg", sep=""), width =5, height = 3, unit ='in', res=1000)
  colors = make_col_scheme(chars = c('A', 'C', 'G', 'T'), cols = c('red2', 'blue4', 'darkorange2', 'green4'))
  g=ggseqlogo(as.character(df[,"plotcol"]), method='bits',seq_type='dna', col_scheme=colors ) + theme_classic() + xlab("") + ylab("") +
	  scale_y_continuous(limits = c(0, 2), expand = c(0.02, 0)) + scale_x_continuous(breaks = seq(1,7,1), labels = seq(-3,3,1), expand = c(0.02, 0)) +
	  theme(
          axis.line = element_line(color = "black", size = 1), axis.ticks = element_line(color = "black", size = 1), axis.ticks.length = unit(0.5, "cm"), axis.text = element_text(color = "black", size = 50), plot.margin = unit(c(0.5,0, 0, 0), "cm")
	  )
  return(g)
}