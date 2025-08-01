#!/bin/bash
#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 4
#SBATCH --mem=4GB
#SBATCH --job-name=TrimMaster
###################################################
## Sanity Check for command line parameters
## Can be skipped but variables will be referenced
## later in the script. checks only if existing
## exits for essential parameters
###################################################################################################
if [ -z ${PFAD+x} ]
then
    echo -e "Pathname parameter to exome folders is not set\nExport the proper var during the sbatch script submission with --export=ALL,PFAD='/path/to/folders' as first argument"
    exit
else
    echo -e "Path to exome files is set to '$PFAD'"
fi
if [ -z ${OUTPUT_DIR+x} ]
then
    echo -e "No output directory specified."
    exit
else
    echo -e "Set output directory to $OUTPUT_DIR"
fi

if [ -z ${FolderList+x} ]
then
    echo -e "No path to the folder list \n should be --export=ALL,FolderList='/path/to/folder_list.txt' as first argument\n"
    exit
else
    echo "Current folder list path is '$FolderList'"
fi
if [ -z ${TRIM_TECH_PARAMETERS+x} ]
then
    echo -e "No parameters set \n Choose correct sequencing thechnology\n"
    exit
else
    echo "Current parameters are: '$TRIM_TECH_PARAMETERS'"
fi
if [ -z ${TECH+x} ]
then
    echo -e "No technology  set \n"
    exit
else
    echo "Current sequencing technology is: '$TECH'"
fi

###################################################################################################
if [ -d $OUTPUT_DIR ]
then
    echo -e "Output directory already present, no further action done. Existing files with same names will be overwritten."
else
    echo -e "Creating folder $OUTPUT_DIR"
    mkdir -p $OUTPUT_DIR
fi


echo -e "Read all folder specified in  $FolderList file\n"
readarray -t LINES < "$FolderList"
for sampleID in "${LINES[@]}"; do
    echo -e "\n Performing Trimming in sampleID =  $sampleID \n"
    #Error and ouput names file with the sampleID name for better debugging
    var1="Trim-e-${sampleID}.txt"
    var2="Trim-o-${sampleID}.txt"
    sbatch --error=$var1 --output=$var2 --export=ALL,PFAD=$PFAD,TECH="$TECH",sampleID=$sampleID,OUTPUT_DIR=$OUTPUT_DIR,TRIM_TECH_PARAMETERS="$TRIM_TECH_PARAMETERS" ./scripts/shell/Trim.sh
done

exit
