#! usr/bin/env bash
salloc -A gts-fstorici3-biocluster -phive -N1 --ntasks-per-node=24 -t1:00:00
source ~/.bash_profile

for file in $(ls *.bed)
do
grep 'chrM' $file >mito/$(basename $file .bed)_mito.bed
grep -v 'chrM' $file > nucl/$(basename $file .bed)_nucl.bed 
done

#Bed files after mismatch analysis and poly N (dN tailing) mismatch removal uptill 4 nucleotides

conda activate r_env
#inputloc=$(pwd)/mito
inputloc=$(pwd)/nucl
#inputloc=$(pwd)
resultsloc=$inputloc/hotspots
#fai=$saccergenomechrM
fai=$saccergenomenucl
#fai=/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38-mito.fa.fai

mkdir $resultsloc $resultsloc/out $resultsloc/out/ppthresh $resultsloc/out/percent
############################### Getting counts / Hotspots #####################################

for thresh in 0.01; do
for file in $(ls $inputloc/FS3*.bed)
do
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/hotspots_thresh.R -b $file -g $fai -r BSgenome.Scerevisiae.UCSC.sacCer2 -t $thresh -o $resultsloc 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/hotspots_thresh.R -b $file -g $fai -r BSgenome.Hsapiens.UCSC.hg38 -t $thresh -o $resultsloc &
done
done
wait

#Getting counts for hotspots
for file in $(ls *percentthresh_0.01*); do awk '{sum+=$8} END {print sum}' $file; done
for file in $(ls *percentthresh_0.01*); do echo $file; done
############################### Common hotspots #####################################
for file in $(ls *files); do
#Rscript /storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/df_matrix.R -f $file -c 8 -s -o ${file}_common_counts.tsv &
Rscript /storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/df_matrix.R -f $file -c 10 -s -o ${file}_common_EF.tsv &
done

thresh=0.49

for file in $(ls *files); do
Rscript /storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/df_matrix.R -f $file -a -t $thresh -c 8 -s -o ${file}_${thresh}_common_counts.tsv &
Rscript /storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/df_matrix.R -f $file -a -t $thresh -c 10 -s -o ${file}_${thresh}_common_EF.tsv &
done

#Properties of common hotspots
wc -l *_common_counts.tsv
for file in $(ls *_common_counts.tsv); do
echo $file
cut -f4 $file | sort | uniq -c
done

wc -l *${thresh}_common_counts.tsv
for file in $(ls *${thresh}_common_counts.tsv); do
echo $file
cut -f4 $file | sort | uniq -c
done

################################ Sequence module  #####################################
#RMres='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_2021/results'
#for file in $(ls FS33*.bed)
#do
#lib=$(echo $file | tr '/' '\n' | tail -1 | tr '.' '\n' | head -1 | tr '_' '\n' | head -1)
#cp $file $RMres/$lib/coordinate30/$lib.bed
#cut -f1,2,3,6 $RMres/$lib/coordinate30/$lib.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' > $RMres/$lib/coordinate30/$lib.counts.tab
#done
################################ Getting MEME for expanded regions maybe! #####################################


################################ Checking TSS regions for yeast #####################################
cd EF


for file in $(ls ../nucl/*.bed); do bash ~/p-fstorici3-0/rich_project_bio-storici/bin/DIRR/scripts/bedtoEF.sh $file $saccergenomenucl 200; done
cd both_strands
rename "_nucl.bw" "" *.bw


#!usr/bin/env bash

#Go to the folder where bw files are
in='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/AGS/subnfilt/EF/both_strands'
out='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/AGS/subnfilt/EF/both_strands/output'
mkdir $out
files='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/AGS/subnfilt/EF/both_strands/files'

#while read p; do file=$(echo $p | awk '{print $1}'); replacement=$(echo $p | awk '{print $3}'); cp $file $replacement; done <$files
#cut -f4 $files | uniq

cd $in
function deeptoolsgene {

groupfiles=$(grep $group $files | cut -f3 | tr '\n' ' ')
bref=$(basename $ref)
base="${bref%.*}"

computeMatrix scale-regions -R $ref --beforeRegionStartLength 1000 --afterRegionStartLength 1000 -S $groupfiles -o $out/${base}_${group}.gz --outFileSortedRegions $out/${base}_${group}.bed --numberOfProcessors max --regionBodyLength 1600 -bs 200 --averageTypeBins sum  --missingDataAsZero 

plotProfile -m $out/${base}_${group}.gz -out $out/${base}_${group}.png --numPlotsPerRow 2 --perGroup --yMin 0  --yMax 4 --legendLocation upper-right --endLabel 'TTS'
}

ref=/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Yeast_all/annotationdb/annotmeta.tsv
for ref in /storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Yeast_all/annotationdb/protein_coding_genes /storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Yeast_all/annotationdb/non_coding_genes
do
for group in BY4742-WT BY4742-rnh201-KO BY4742-rnh201-G42S BY4742-rnh203-K46W
do
deeptoolsgene
done
done

