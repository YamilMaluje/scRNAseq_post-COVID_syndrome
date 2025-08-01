#! /bin/bash
#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 18
#SBATCH --mem=8GB
#SBATCH --job-name=QCFastq

###################################################################################################
# Asses the quality of the FASTQ reads with FASTQC
module load singularity-sylabs/v3.8.1
#Activate the singularity shell
singularity shell /data/sb_service_01/YM_containers/scRNAseq_v1.4.sif


###################################################################################################
if [ -z ${STEP+x} ]
then
    echo -e "NO Analysis step assigned closing pipeline..."
    exit
else 
    echo "Previous step: $STEP "
fi

if [ ${SATUS} -eq 0 ]
then
    echo -e "Error in the Previous step: $STEP"
    exit
elif [ ${SATUS} -eq 1 ]
    echo -e "OK in the Previous step: $STEP"
else 
    echo "STATUS NOT SET = $SATUS"
fi

if [ -z ${ORIG_PATH+x} ]
then
    echo -e "NO project path set..."
    exit
else 
    echo "Current project path is '$ORIG_PATH'"
fi

###################################################################################################
#Obtaining project path
ORIG_PATH=$(pwd)
echo -e "Project path: $ORIG_PATH"

#FASTQC
###################################################################################################
# Check if 00_Integrity_check_files.txt file exists
REPORT_FILE_2=( $( find ${ORIG_PATH}/00_Quality_Check -type f  \( -iname 00_Integrity_check_files.txt \) ))

# Create fastqc file output if 00_Integrity_check_files.txt file exists
if [[ ! -f $REPORT_FILE_2 ]]; then
  echo "File $REPORT_FILE_2 not found!"
  exit 1
else
  dirname_path=( $(dirname ${REPORT_FILE_2} ))
  #Create a folder for the report output
  mkdir -p ${dirname_path}/01_fastqc
fi
#Find all fastq files
FASTQ_FILES=$(  find ${ORIG_PATH}/00_fastq -type f  \( -name *.fastq -o -name *.fastq.gz -o -name *.fq -o -name *.fq.gz  \) )
#Running fastqc on the found fastq files
fastqc -t 18 -o ${dirname_path}/01_fastqc ${FASTQ_FILES}
exit




#Multiqc   
###################################################################################################
#Run MULTIQC on the output of FASTQC 

singularity shell /data/sb_service_01/YM_containers/multiqc-1.9.sif
#Obtaining project path
ORIG_PATH=$(pwd)

FASTQC_FILES=$(  find ${ORIG_PATH}/00_Quality_Check/01_fastqc -type f  \( -name *R1*_fastqc.zip -o -name *R2*_fastqc.zip -o -name *_1_fastqc.zip -o -name *_2_fastqc.zip \) )

# Create 02_multiqc_fastq  if fastqc worked
if [[ -z "$FASTQC_FILES" ]]; then
  echo "FASTQC did not work"
  exit 1
else
  echo "File $FASTQC_FILES found!"
  #Create a folder for the report output
  OUT_DIR=${ORIG_PATH}/00_Quality_Check/02_multiqc_fastq
  mkdir -p $OUT_DIR
fi

#running multiqc
if [[ -n ${OUT_DIR} ]]; then
  multiqc -o $OUT_DIR $FASTQC_FILES
else
  exit
fi

#check if the multiqc file was created
MULTIQC_FILES=$(  find ${OUT_DIR} -type f  \( -name multiqc_report.html \) )
if [[ -n ${MULTIQC_FILES} ]]; then
  echo "Multiqc runned well"
else
  exit
fi
#exit singularity image
exit
exit
