# rNMP_hotspots

Finding common and highly incorporated rNMP locations using rNMP Enrichment Factor

## Installation

Use [conda](https://pip.pypa.io/en/stable/](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to create enviroment from r_env.yml file.

```bash
conda env create -n r_env --file r_env.yml
```

## Usage

```bash
# Activate Enviroment
conda activate r_env

# Define Variables
inputloc=$(pwd)/mito #Location of bed files
resultsloc=$inputloc/hotspots
scripts=rNMP_hotpots/scripts #location of this repository
fai=~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38-mito.fa.fai #full path of reference genome .fai file
thresh=0.01 #Fraction of locaitons with highest rNMP counts to be filtered


# One script run for getting rNMP counts, enrichment Factor, filter hotspots from each bed file and generate ggseqlogo plots for hotspots (filtered based on threshold provided)

mkdir $resultsloc $resultsloc/out $resultsloc/out/ppthresh $resultsloc/out/percent
for file in $(ls $inputloc/FS3*.bed)
do
Rscript $scripts/hotspots_thresh.R -b $file -g $fai -r BSgenome.Hsapiens.UCSC.hg38 -t $thresh -o $resultsloc &
done
wait


# Get counts and Enrichment Factor matrices based on cell type 
mat=$scripts/matrix
Rscript $scripts/df_matrix.R -f $mat -a -t $thresh -c 10 -s -o ${file}_${thresh}_common_EF.tsv
Rscript $scripts/df_matrix.R -f $mat -a -t $thresh -c 8 -s -o ${file}_${thresh}_common_counts.tsv

```


## License

[MIT](https://github.com/DKundnani/rNMP_hotspots/blob/main/LICENSE)
