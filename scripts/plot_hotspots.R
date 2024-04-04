#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", args[grep(file.arg.name, args)])
script.basename <- dirname(script.name)
other.name <- file.path(script.basename, "hotspot_analysis.R")
#other.name <- file.path("/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR", "hotspot_analysis.R")
source(other.name)

#Taking in arguments
library("optparse")
option_list = list(
  make_option(c("-m", "--mat"), type="character", default=NULL, help="ribo bed file with col: chr start stop . . strand",metavar="character"),
  make_option(c("-c", "--col"), action = "store_true", default = FALSE, help="if the matrix has header or not",metavar="character"),
  make_option(c("-g", "--gen"), type="character", default=NULL, help="genome (.fai) file", metavar="character"),
  make_option(c("-t", "--thresh"), type="character", default=NULL, help="threshold to select, if ", metavar="character"),
  make_option(c("-r", "--ref"), type="character", default='BSgenome.Hsapiens.UCSC.hg38', help="BS genome name to be used for reference", metavar="integer"),
  make_option(c("-v", "--vis"), action = "store_true", default = FALSE, help="visualize or not",metavar="character"),
  make_option(c("-o", "--output_prefix"), type="character", default="out", help="output file name [default %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#Checking for required input
writeLines("\n...Checking input...\n")
if (is.null(opt$mat)){
  print_help(opt_parser)
  stop("Please specify matrix file.n", call.=FALSE)
}
if (is.null(opt$gen)){
  print_help(opt_parser)
  stop("Please specify genome (.fai) file .n", call.=FALSE)
}

mat=opt$mat; writeLines(mat)
c=opt$col; print(c) #cov='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/AGS/DNAseq/aligned/Y1.cov'
gen=opt$gen; writeLines(gen) #gen='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/sacCer2/sacCer2-nucl.fa.fai'
out=opt$output_prefix; writeLines(out)
thresh=as.numeric(opt$thresh)
ref=opt$ref
vis=opt$vis
#####Main Script
#1.Calling required Packages
writeLines("\n...Calling required Package(s)...\n")
suppressMessages(library(tools,quietly = T)) #for getting file basename
suppressMessages(library(bedr, quiet = T))
suppressMessages(library(regioneR,quiet = T))
suppressMessages(library(dplyr,quiet = T))
suppressMessages(library(data.table,quiet = T))
suppressMessages(library(stringr,quiet = T))
suppressMessages(library(zoo,quiet = T))
suppressMessages(library(karyoploteR, quiet = T))
suppressMessages(library(httpgd, quiet = T))
pre=file_path_sans_ext(basename(mat))

#2.Defining functions if any
#writeLines("\n...Defining required functions...\n")

#3. Preprocessing input files
#writeLines("\n...Processing Input files...\n")
genome<-read.table(gen,sep="\t",header=F)
genome=genome[genome[1]!="2micron",]
custom.genome <- toGRanges(data.frame(chr=genome[1], start=rep(1, nrow(genome)), end=genome[2]))
ribo=c("A", "C", "G", "U")
col=c("#D55E00","#0072B2", "#DFAC00","#009E73") #rA,rC, rG, rU
#Reading matrix
if (c) {
  stats<-read.table(mat,sep="\t",header=T); stats=unique(stats)
} else {
  stats<-read.table(mat,sep="\t",header=F); stats=unique(stats)
}

#4.Main code
#writeLines("\n...Executing Main code...\n")
n=ncol(stats)
#samt=round((n-8)*0.5)
samt=2
colnames(stats)[1:6]<-c("chr", "start", "stop", "V4", "V5", "strand")
stats$percrank<-(ecdf(stats$average)(stats$average)+ecdf(stats$samplexists)(stats$samplexists))/2
filtstats=stats[stats[n-1]>=samt,]
if (thresh<1) {
  filtstats=filtstats[filtstats$percrank>=(1-thresh),]
} else {
  filtstats=head(filtstats[order(filtstats$percrank, decreasing=TRUE),],thresh)
}

#filtstats <- stats[order(stats[,(n-1)],stats[,n], na.last = TRUE, decreasing = TRUE ),][1:top,]
final=get_flank_seq(filtstats,3,3, ref)

writeLines(mat)
print(table(final$ribo))

final$col<-rep("", nrow(final))
for (r in seq(1,4)) {
  if (nrow(final[final$ribo==ribo[r],])) {
    final[final$ribo==ribo[r],]$col<-col[r]
  }
}

#final$seqU<-paste(substr(unlist(final$seqU),1,3), rep("r", nrow(final)),substr(unlist(final$seqU),4,4), substr(unlist(final$seqU),5,7), sep="")
final.df<-final[,c("chr","start","stop", "ribo","seqU","strand",colnames(final)[n-1],colnames(final)[n],"seq"),]
write.table(final.df, paste(str_split(mat, pattern="[.]")[[1]][1],"_",thresh,"_top.tsv", sep=""), sep='\t',quote = F, col.names=F, row.names=F)


if (vis) {
  statspos<-final[final$strand == '+',]
  statsneg<-final[final$strand == '-',]
  statspos.gr<-makeGRangesFromDataFrame(statspos, keep.extra.columns=T, starts.in.df.are.0based=T)
  statsneg.gr<-makeGRangesFromDataFrame(statsneg, keep.extra.columns=T, starts.in.df.are.0based=T)
  
  if (nrow(genome)>1) {
    #jpeg(paste(str_split(mat, pattern="[.]")[[1]][1],"_vis_label.jpeg", sep=""), width=25, height = 30, unit ='in', res=1000)
    #pp <- getDefaultPlotParams(plot.type=2)
    #pp$ideogramheight=8
    #pp$data1height=500
    #pp$data2height=500
    #kp <- plotKaryotype(genome = custom.genome, plot.type=2, plot.params = pp) #levels(seqnames(regions))) # nolint
    #kpAddCytobandsAsLine(kp)
    #kpAddBaseNumbers(kp, tick.dist = 5e5, cex=1, minor.tick.dist = 1e5)
    #kpPoints(kp, data=statspos.gr, y=statspos.gr$average*max(stats$average)/100, cex=2, data.panel = 1, ymax=1, ymin=0, col=statspos.gr$col)
    #kpPlotMarkers(kp, data=statspos.gr, labels = paste(statspos$stop, statspos$ribo, sep=":"), font = 2, text.orientation = "vertical", data.panel = 1, marker.parts = c(0.3, 0.4, 0.3), label.color = statspos.gr$col, adjust.label.position=TRUE, label.dist = 0.002, y=0.4, max.iter=1000)
    #kpPlotMarkers(kp, data=statsneg.gr, labels = paste(statsneg$stop, statsneg$ribo, sep=":"), font = 2, text.orientation = "vertical", data.panel = 2, marker.parts = c(0.3, 0.4, 0.4), label.color = statsneg.gr$col,  adjust.label.position=TRUE, label.dist = 0.002, y=0.4, max.iter=1000)
    #kpPoints(kp, data=statspos.gr, y=0.15, cex=1.5, data.panel = 1, ymax=1, ymin=0, col=statspos.gr$col)
    #dev.off()

    #jpeg(paste(str_split(mat, pattern="[.]")[[1]][1],"_vis_label_hori.jpeg", sep=""), width=10, height = 10, unit ='in', res=300)
    #pp <- getDefaultPlotParams(plot.type=2)
    #pp$ideogramheight=8
    #pp$data1height=350
    #pp$data2height=350
    #kp <- plotKaryotype(genome = custom.genome, plot.type=2, plot.params = pp) #levels(seqnames(regions))) # nolint
    #kpPlotMarkers(kp, data=statspos.gr, labels = paste(statspos$stop, statspos$ribo, sep=":"), font = 2, text.orientation = "horizontal", data.panel = 1, marker.parts = c(0.3, 0.4, 0.3), label.color = statspos.gr$col, adjust.label.position=TRUE, label.dist = 0.002, y=0.4, max.iter=1000)
    #kpPlotMarkers(kp, data=statsneg.gr, labels = paste(statsneg$stop, statsneg$ribo, sep=":"), font = 2, text.orientation = "horizontal", data.panel = 2, marker.parts = c(0.3, 0.4, 0.4), label.color = statsneg.gr$col,  adjust.label.position=TRUE, label.dist = 0.002, y=0.4, max.iter=1000)
    #dev.off()

    jpeg(paste(str_split(mat, pattern="[.]")[[1]][1],"_vis_mark.jpeg", sep=""), width=6, height = 6, unit ='in', res=300)
    pp <- getDefaultPlotParams(plot.type=2)
    pp$ideogramheight=5
    pp$data1height=400
    pp$data2height=400
    kp <- plotKaryotype(genome = custom.genome, plot.type=2, plot.params = pp, labels.plotter = NULL) #levels(seqnames(regions))) # nolint
    kpAddChromosomeNames(kp, cex = 1.2, xoffset=-0.002)
    kpPlotMarkers(kp, data=statspos.gr, labels = statspos$ribo, font = 2, text.orientation = "horizontal", data.panel = 1, label.color = statspos.gr$col,marker.parts =  c(0.3,0.4,0.3),  adjust.label.position=TRUE, label.dist = 0.002, y=0.35, max.iter=1000)
    kpPlotMarkers(kp, data=statsneg.gr, labels = statsneg$ribo , font = 2, text.orientation = "horizontal", data.panel = 2, label.color = statsneg.gr$col,marker.parts =  c(0.3,0.4,0.3),  adjust.label.position=TRUE,label.dist = 0.002, y=0.35, max.iter=1000)
    dev.off()

  } else {
    #jpeg(paste(str_split(mat, pattern="[.]")[[1]][1],"_vis_label.jpeg", sep=""), width=5, height = 3, unit ='in', res=300)
    #pp <- getDefaultPlotParams(plot.type=2)
    #pp$ideogramheight=2
    #pp$data1height=500
    #pp$data2height=500
    ##kp <- plotKaryotype(genome = custom.genome, plot.type=2, plot.params = pp) #levels(seqnames(regions))) # nolint
    #kpPlotMarkers(kp, data=statspos.gr, labels = paste(statspos$stop, statspos$ribo, sep=":"), font = 2,text.orientation = "vertical", data.panel = 1, marker.parts = c(0.2, 0.6, 0.2), label.color = statspos.gr$col,  label.dist = 0.005, y=0.1, clipping = TRUE)
    #kpPlotMarkers(kp, data=statsneg.gr, labels = paste(statsneg$stop, statsneg$ribo, sep=":"), font = 2,text.orientation = "vertical", data.panel = 2, marker.parts = c(0.2, 0.6, 0.2), label.color = statsneg.gr$col,  label.dist = 0.005, y=0.1, clipping = TRUE)
    #dev.off()
    
    jpeg(paste(str_split(mat, pattern="[.]")[[1]][1],"_vis_mark.jpeg", sep=""), width=5, height = 3, unit ='in', res=300)
    pp <- getDefaultPlotParams(plot.type=2)
    pp$ideogramheight=2
    pp$data1height=500
    pp$data2height=500
    kp <- plotKaryotype(genome = custom.genome, plot.type=2, plot.params = pp) 
    kpPlotMarkers(kp, data=statspos.gr, labels = statspos$ribo, font = 2,text.orientation = "horizontal", data.panel = 1, label.color = statspos.gr$col,marker.parts =  c(0.6,0.2,0.2),  adjust.label.position=TRUE, label.dist = 0.005, y=0.2, max.iter=1000)
    kpPlotMarkers(kp, data=statsneg.gr, labels = statsneg$ribo, font = 2,text.orientation = "horizontal", data.panel = 2, label.color = statsneg.gr$col,marker.parts =  c(0.6,0.2,0.2),  adjust.label.position=TRUE,label.dist = 0.005, y=0.2, max.iter=1000)
    dev.off()
  }
}

#writeLines("\nJob Finished.\n")
#writeLines(paste("Output folder name:", out))
#writeLines(paste("Output columns:", c("chr", "start", "loc", "counts", "round_norm_counts", "strand", "norm_count", "EF")))
#write.table(final.df, file=paste(out,"/",pre,".norm.counts", sep=""), sep = "\t",quote=F,row.names=F, col.names=F)
#write.table(final.df[final.df$round_norm_counts>0,], file=paste(out,"/",pre,".bed", sep=""), sep = "\t",quote=F,row.names=F, col.names=F)
