#! usr/bin/bash

get_input () {
	
	while getopts "e:i:fchisdb" option
	do 
		case ${option} in
			e) env=${OPTARG}; ;;
			o) resultsloc=${OPTARG}; ;;
			r) ref=${OPTARG}; ;;
			t) thresh=${OPTARG}; ;;
      b) bed=${OPTARG}; ;;
			h) echo "usage: bash lib_hotspots.sh -env r_env -i  -c -o output -b bedfile(s)
					[-e <enviroment to be activated>] 
					[-o <output folder>] 
					[-r <reference fai file>] 
					[-t <threshold in fraction(percentile) of unique locations to be considered as hotspots>]
					[-b <annotated files will be separated based on nucleus and mitochondria>]
					[-h] help option, gives usage information 
				";;
			*) 
	
		esac
	done

	
	#flag check
	if [[ $reference == "" ]]; then
	echo "Please provide a reference"
	exit 1
	fi

	if [[ $b == "" ]]; then
	echo "Please provide bed file "
	exit 1
	fi

	if [[ $resultsloc == "" ]] ; then
	echo "Using current directory as output directory"
	fi
  
  if [[ $env == "" ]] ; then
	echo "Assuming the enviroment is already activated or required dependencies are readily availbel for the script to use"
  else
  conda activate $env
  echo "Enviroment activated"
	fi


}


initiate () {
	#check for reference annotation and bed files 
	if [[ -f "$reference" ]]; then echo "$reference found"; else echo "annotation file not found in the specified location"; fi
	echo 
	#for file in "$files"
	for file in ${bed} 
	do
	if [[ -f $file ]];then echo "$file found"; else "Bed files not found in the current locations"; fi 
	done
	echo 
	#create output directory
	if [[ -d $output ]]; then echo "Removing already existing output directory"; rm -r $output;mkdir $output; [[ -d "annotations" ]] && echo "New output directory created"  ; else mkdir $output; [[ -d $output ]] && echo "output directory created"; fi
	echo
	#check if the reference annotation format is compatible for strandness
	if [[ "$strand" == 1 ]]; then col6=$(cut -f6 "$reference" | head -1); if [[ ${col6} =~ "+" ]] || [[ "${col6}" =~ "-" ]]; then echo "Reference annotation format looks good"; else echo "Reference format not compatible. Please have bed6 file with strand information in the 6th column"; fi; else echo "No strand information needed"; fi
  	if [[ $divide == 1 ]]; then echo "Diving by organelle"; else echo "Not diving by organelle"; fi
}

gethotspots () {
thresh=${1} #0.01
mkdir $resultsloc $resultsloc/out $resultsloc/out/ppthresh $resultsloc/out/percent
for file in $(ls $inputloc/*.bed)
do
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/hotspots_thresh.R -b $file -g $fai -t $thresh -o $resultsloc &
done
wait

}
############################### Getting counts / Hotspots #####################################





thresh=2
mkdir $resultsloc $resultsloc/out $resultsloc/out/ppthresh $resultsloc/out/percent
for file in $(ls $inputloc/*.bed)
do
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/hotspots_thresh.R -b $file -g $fai -t $thresh -o $resultsloc &
done
wait


main() {
	get_input "$@"
	initiate 
	annotating 
  split
  allcounts
}

# Calling the main function
main "$@"






