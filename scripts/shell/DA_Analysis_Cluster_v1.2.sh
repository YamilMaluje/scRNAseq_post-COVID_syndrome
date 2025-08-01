#!/bin/bash
#SBATCH --partition=shortterm
#SBATCH --job-name=DA_scRNAseq
#SBATCH --no-requeue
#SBATCH --mem=48GB
#SBATCH -c 12
###
###################################################################################################
start=`date +%s`
## Load singularity module
#module load singularity-sylabs/v3.8.1
module load singularity/v3.8.1
echo -e "Differentially Abundance script Initializing... \n\n\n "
###################################################################################################
if [ -z ${SEU_OBJ_PATH+x} ]
then
    echo -e "No seurat object path found..."
    exit
else
    echo -e "Seurat object path is:  '$SEU_OBJ_PATH'"
fi
if [ -z ${OUTPUT_DIR+x} ]
then
    echo -e "No Output path provided...."
    exit
else
    echo "Output path is :  '$OUTPUT_DIR'"
fi
if [ -z ${VARIABLE_COMPARISON+x} ]
then
    echo -e "No Variable of comparison provided provided...."
    echo -e "If the user wants to compare the Differencially abundance of Cancer cells or not the following promt must be created \n"
    echo -e "Cancer:0:1"
    echo -e "This will split the seurat object by the metadata 'Cancer' and will dispay the cells that are tagged with '0' as Blue and tagged with '1' as Red"
    echo -e "Note: The promt must be separated by ':' with no spaces"
    exit
else
    echo "Variable of comparison is :  '$VARIABLE_COMPARISON'"
fi
if [ -z ${REDUCTIONS+x} ]
then
    echo -e "No Variable of comparison provided provided...."
    echo -e "Dimensional reduction to use for plots (i.e. UMAP, TSNE), could be a string or a string vector\n"
    exit
else
    echo "Dimensional Reduction for plots :  '$REDUCTIONS'"
fi


###################################################################################################
echo "Running R script "
singularity exec /data/sb_service_01/YM_containers/scRNAseq_v1.4.sif Rscript /data/sb_service_01/YM_containers/scripts/R/DA_Analysis_Cluster_v1.2.R $SEU_OBJ_PATH $OUTPUT_DIR $VARIABLE_COMPARISON $REDUCTIONS
end=`date +%s`
echo -e "Everything Done ..."
echo -e "Execution time was `expr $end - $start` seconds."
exit
