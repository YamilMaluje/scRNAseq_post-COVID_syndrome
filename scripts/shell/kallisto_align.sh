#! /bin/bash
#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 24
#SBATCH --mem=96GB
#SBATCH --job-name=KallAlign

########################################################
## Single Sample Kallisto Single Cell Alignment
## Start multiple instances if you have several samples


###################################################
## Load singularity module
module load singularity-sylabs/v3.8.1
echo -e "Kallisto Align script Initializing... \n\n\n "

mkdir $SCRATCH/index
mkdir $SCRATCH/index/inp
mkdir $SCRATCH/index/assembled
mkdir $SCRATCH/source/
mkdir $SCRATCH/alignment/
mkdir $OutputDir/$sampleID/


ORIG_PATH=$(pwd)
###################################################
## Sanity Check for command line parameters
## Can be skipped but variables will be referenced
## later in the script. checks only if existing
## exits for essential parameters


if [ -z ${OUTPUT_DIR+x} ]
then
	echo -e "No path to output provided\n_Export the proper var during the sbatch script submission with --export=ALL,OutputDir='/path/to/Outputdir/' as first argument\n__Escaping script"
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
	INDEX=${Index_Assemb}
fi

if [ -z ${T2g+x} ]
then
	echo -e "Transcript to gene file is missing \n If not provided, build it with option 'Index generator'"
	exit
else 
	echo -e "Path to kallisto transcript to gene '$T2g'"
	T2G=${T2g}
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

if [ -z ${sampleID+x} ]
then
    echo -e "No sample added \n"
    exit
else
    echo "Current sample: '$sampleID'"
fi


#################################################
if [ -d $OUTPUT_DIR/$sampleID ]
then
    echo -e "Output directory already present, no further action done. Existing files with same names will be overwritten."
else
	echo -e "Creating folder $OUTPUT_DIR/$sampleID"
    mkdir -p $OUTPUT_DIR/$sampleID/"cells_x_genes"
fi

rsync -av $Whitelist $SCRATCH/index/whitelist.txt


echo -e "List of technologies suported:"
singularity exec /data/sb_service_01/YM_containers/scRNAseq_v1.4.sif kallisto bus -l

case $TECH in 
	"10x Genomics v2")
		echo -e "10xv2 technology was set "
		AUX="10xv2"
	;;
	"10x Genomics v3")
		AUX="10xv3"
		echo -e "10xv3 technology was set "
	;;
	"Bd Rhapsody")
		AUX="0,0,9,0,21,30,0,43,52:0,52,60:1,0,0"
		echo -e "Bd Rhapsody technology was set "
	;;
	"Cscope")
		AUX="0,0,8,0,24,32,0,48,56:0,57,69:1,0,0"
		echo -e "Cscope technology was set "
	;;
	"scopeV3.0.1")
		AUX="0,0,9,0,25,34,0,50,59:0,60,72:1,0,0"
		echo -e "Cscope V3.0.1 technology was set"
		echo -e "C9 L16 C9 L16 C9 L1 U12 T18"
	;;
	*)
		echo -e "Other Technology is not supported."
	;;
esac




## Reference Genome
##INDEX=( $(  find ${Index_Assemb} -type f  \( -iname *.idx \) ))

## Reference Genome
## If no Reference Genome supplied
## one from Repository will be used
rsync -av $INDEX $SCRATCH/index/assembled/transcriptome.idx

## Transcript to Gene textfile
##T2G=( $(  find ${Index_Assemb} -type f  \( -iname transcripts_to_genes.txt \) ))
## Transcript to Gene textfile
rsync -av $T2G $SCRATCH/index/t2g.txt

#grab the fastq files path
R1=( $(  find ${ReadDir} -type f  \( -iname *${sampleID}*1.fastq -o -iname *${sampleID}*1.fastq.gz -o -iname *${sampleID}*1.fq -o -iname *${sampleID}*1.fq.gz  \) ))
R2=( $(  find ${ReadDir} -type f  \( -iname *${sampleID}*2.fastq -o -iname *${sampleID}*2.fastq.gz -o -iname *${sampleID}*2.fq -o -iname *${sampleID}*2.fq.gz  \) ))

echo -e "Index path: ${INDEX}"
echo -e "Transcript to gene path: ${T2G}"
echo -e "Read1: ${R1}"
echo -e "Read2: ${R2}"


## Read in Sample from path and unzip
cp -a $R1 $SCRATCH/source/
cp -a $R2 $SCRATCH/source/
echo "This are all the files coppied: "
ls $SCRATCH/source
gunzip $SCRATCH/source/*.gz
echo "This are all the files coppied UNZIPPED: "
ls $SCRATCH/source

## Set path variables for R1 and R2
READ1=( $(  find $SCRATCH/source/ -type f  \( -iname *${sampleID}*1.fastq -o -iname *${sampleID}*1.fq  \) ))
READ2=( $(  find $SCRATCH/source/ -type f  \( -iname *${sampleID}*2.fastq -o -iname *${sampleID}*2.fq  \) ))
echo "Read one containing the barcodes in path $READ1"
echo "Read two containing the transcript in path $READ2"


#cd $SCRATCH/alignment/

#RUNNUING KALLISTO BUS
COMMAND=""
COMMAND+="kallisto bus "
COMMAND+="-i $SCRATCH/index/assembled/transcriptome.idx " 
COMMAND+="-o $SCRATCH/alignment/ "
COMMAND+="-x ${AUX} "
COMMAND+="-t 24 "
COMMAND+="${READ1} "
COMMAND+="${READ2}"

echo -e "\n"
echo "This is the command for Kallisto BUS: $COMMAND"
echo -e "\n"

echo "Starting Kallisto Bus..."
singularity exec /data/sb_service_01/YM_containers/scRNAseq_v1.4.sif $COMMAND
echo "Finished Kallisto Bus..."


mkdir $SCRATCH/alignment/genecount/ $SCRATCH/alignment/tmp/
## create matrices and final output based on the busfile
COMMAND_CORRECT=""
COMMAND_CORRECT+="bustools correct "
COMMAND_CORRECT+="-w $SCRATCH/index/whitelist.txt "
COMMAND_CORRECT+="-o $SCRATCH/alignment/temp_output.bus "
COMMAND_CORRECT+=" $SCRATCH/alignment/output.bus"

singularity exec /data/sb_service_01/YM_containers/scRNAseq_v1.4.sif $COMMAND_CORRECT 


COMMAND_SORT=""
COMMAND_SORT+="bustools sort "
COMMAND_SORT+="-T $SCRATCH/alignment/tmp/ "
COMMAND_SORT+="-t 4 "
COMMAND_SORT+="-o $SCRATCH/alignment/sort_output.bus "
COMMAND_SORT+=" $SCRATCH/alignment/temp_output.bus "

singularity exec /data/sb_service_01/YM_containers/scRNAseq_v1.4.sif $COMMAND_SORT 


COMMAND_COUNT=""
COMMAND_COUNT+="bustools count "
COMMAND_COUNT+="-o $SCRATCH/alignment/genecount/cells_x_genes "
COMMAND_COUNT+="-g $SCRATCH/index/t2g.txt "
COMMAND_COUNT+="-e $SCRATCH/alignment/matrix.ec "
COMMAND_COUNT+="-t $SCRATCH/alignment/transcripts.txt "
COMMAND_COUNT+="--genecounts "
COMMAND_COUNT+=" $SCRATCH/alignment/sort_output.bus "

singularity exec /data/sb_service_01/YM_containers/scRNAseq_v1.4.sif  $COMMAND_COUNT

echo "Starting bustools correct | bustools sort | bustools count..."

### move output to outputDir
cp $SCRATCH/alignment/genecount/cells_x_genes.* $OUTPUT_DIR/$sampleID/"cells_x_genes"


mv "KallAlign-"*"-${sampleID}.txt" ${OUTPUT_DIR}
rm -rf $SCRATCH/*

echo -e "Everything Done correctly...\n"
exit
