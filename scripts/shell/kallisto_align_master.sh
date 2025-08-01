#! /bin/bash
#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 4
#SBATCH --mem=4GB
#SBATCH --job-name=KallAlignMaster

###################################################
## Sanity Check for command line parameters
## Can be skipped but variables will be referenced
## later in the script. checks only if existing
## exits for essential parameters
if [ -z ${FolderList+x} ]
then
	echo -e "No path to the folder list \n should be --export=ALL,FolderList='/path/to/folder_list.txt/' as first argument\n"
	exit
else 
	echo "Current folder list path is '$FolderList'"
fi

if [ -z ${OUTPUT_DIR+x} ]
then
	echo -e "No path to output provided"
	exit
else 
	echo "Current output path is '$OUTPUT_DIR'"
fi

if [ -z ${Index_Assemb+x} ]
then
	echo -e "Index path or assembly path is missing \n If not provided, build it with option 'Index generator'"
	exit
else 
	echo -e "Path to kallisto index and assembly is set to '$Index_Assemb'"
fi

if [ -z ${T2g+x} ]
then
	echo -e "Transcript to gene file is missing \n If not provided, build it with option 'Index generator'"
	exit
else 
	echo -e "Path to kallisto transcript to gene '$T2g'"
fi

if [ -z ${ReadDir+x} ]
then
	echo -e "No path to input file provided\n_Export the proper var during the sbatch script submission with --export=ALL,ReadDir='/path/to/SampleDir/' as first argument\n__Escaping script"
	exit
else 
	echo "Path to input files is set to '$ReadDir'"
fi

if [ -z ${Whitelist+x} ]
then
	echo -e "Whitelist parameter is not set\n_Export the proper var during the sbatch script submission with --export=ALL,Whitelist='/path/to/whitelist.txt' as first argument"
	exit
else 
	echo -e "Path to whitelist is set to '$Whitelist'"
fi

if [ -z ${TECH+x} ]
then
    echo -e "No technology  set \n"
    exit
else
    echo "Current sequencing technology is: '$TECH'"
fi



echo -e "Read all folder specified in  $FolderList file\n"
readarray -t LINES < "$FolderList"

for folder in "${LINES[@]}"; do 

	echo -e "\n Performing Kallisto alignment on folder: $folder \n"
		
	#Error and ouput names file with the folder name for better debugging 
	var1="KallAlign-e-$folder.txt"
	var2="KallAlign-o-$folder.txt"
	
	sbatch --error=$var1 --output=$var2 --export=ALL,TECH="$TECH",Whitelist=$Whitelist,ReadDir=$ReadDir,sampleID=$folder,OUTPUT_DIR=$OUTPUT_DIR,Index_Assemb=$Index_Assemb,T2g=$T2g ./scripts/shell/kallisto_align.sh
	
done


exit

