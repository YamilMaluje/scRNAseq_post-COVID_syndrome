#!/bin/bash
#SBATCH --partition=shortterm
#SBATCH --job-name=Trim
#SBATCH --no-requeue
#SBATCH --mem=48GB
#SBATCH -c 12
###
###################################################################################################
## Load singularity module
module load singularity-sylabs/v3.8.1
echo -e "Trimming script Initializing... \n\n\n "


mkdir $SCRATCH/inp
mkdir $SCRATCH/out

###################################################################################################
if [ -z ${PFAD+x} ]
then
echo -e "Pathname parameter to exome folders is not set\nExport the proper var during the sbatch script submission with --export=ALL,PFAD='/path/to/folders' as first argument"
exit
else
echo -e "Path to exome files is set to '$PFAD'"
fi
if [ -z ${sampleID+x} ]
then
    echo -e "Providing a sample identifier is mandatory\n_Provide sample identifier in the form of >s1< for the first sample etc\n__Export the proper var during the sbatch script submission with --export=ALL,SampleID='sampleID' as first argument\n___Escaping script"
    exit
else
    echo "Current sample processed is '$sampleID'"
fi
if [ -z ${OUTPUT_DIR+x} ]
then
    echo -e "No output directory specified."
    exit
else
    echo -e "Set output directory to $OUTPUT_DIR"
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
if [ -d $OUTPUT_DIR/$sampleID ]
then
    echo -e "Output directory already present, no further action done. Existing files with same names will be overwritten."
else
    echo -e "Creating folder $OUTPUT_DIR"
    mkdir -p $OUTPUT_DIR/$sampleID
fi



case $TECH in 
    "10xv3")
        echo -e "10xv3 technology was set "
        R1=( $(  find ${PFAD} -type f  \( -iname *${sampleID}*R1*.fastq -o -iname *${sampleID}*R1*.fastq.gz -o -iname *${sampleID}*R1*.fq -o -iname *${sampleID}*R1*.fq.gz \) ))
        R2=( $(  find ${PFAD} -type f  \( -iname *${sampleID}*R2*.fastq -o -iname *${sampleID}*R2*.fastq.gz -o -iname *${sampleID}*R2*.fq -o -iname *${sampleID}*R2*.fq.gz \) ))
    ;;
    *)
        echo -e "Other Technology is not supported."
        R1=( $(  find ${PFAD} -type f  \( -iname *${sampleID}*1.fastq -o -iname *${sampleID}*1.fastq.gz -o -iname *${sampleID}*1.fq -o -iname *${sampleID}*1.fq.gz  \) ))
        R2=( $(  find ${PFAD} -type f  \( -iname *${sampleID}*2.fastq -o -iname *${sampleID}*2.fastq.gz -o -iname *${sampleID}*2.fq -o -iname *${sampleID}*2.fq.gz  \) ))
    ;;
esac


echo -e "$sampleID Read 1 =\n ${R1} \n"
echo -e "$sampleID Read 2 =\n ${R2} \n"
temp1=$( basename $R1 )
temp2=$( basename $R2 )
rsync -av $R1 $SCRATCH/inp
rsync -av $R2 $SCRATCH/inp

echo -e "Performing trimming of  R1= ${temp1} file"
echo -e "Performing trimming of  R2= ${temp2} file \n"
COMMAND=""
COMMAND+="/home/fastp  "
COMMAND+="--in1 $SCRATCH/inp/${temp1} "
COMMAND+="--in2 $SCRATCH/inp/${temp2} "
COMMAND+="--out1 $SCRATCH/out/trim.${temp1} "
COMMAND+="--out2 $SCRATCH/out/trim.${temp2} "
COMMAND+="--html $SCRATCH/out/${sampleID}.html "
COMMAND+="$TRIM_TECH_PARAMETERS"
echo -e "\n"
echo "This is the command for fastp: $COMMAND"
echo -e "\n"

echo "Running fastp..."
singularity exec /data/sb_service_01/YM_containers/scRNAseq_v1.4.sif $COMMAND
echo "Finished fastp..."
mv "Trim-"*"-${sampleID}.txt" ${OUTPUT_DIR}
rsync -av $SCRATCH/out/ $OUTPUT_DIR/$sampleID
rm -rf $SCRATCH/*
exit
