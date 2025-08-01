library(Seurat)
library(dplyr)
library(Matrix)
library(DropletUtils)
library(tidyverse)
library(ggpointdensity)
library(scico)
library(scales)
library(BUSpaRse)
library(ggplot2)
library(ggrepel)
library("readxl")
library(patchwork)
################# SETTING VARIABLES ##############################################
orig <- getwd()
#Adjust the single cell path and project path
#proj.path <- paste(orig,"scRNAseq Analysis", "Köhl", "sc_Köhl_mouse", sep = "/")
srcRaw <- paste(orig,"05_DownstreamAnalysis","SelectedGenes", "09.05.2025",sep = "/")
today <- Sys.Date() %>% format(., format="%d.%m.%Y")
output.path <- paste(orig,"05_DownstreamAnalysis","DGE_scRNAseq_Condition",today,sep = "/")
output.path.plots <- paste(output.path, "Plots",sep = "/")
dir.create(output.path.plots,recursive = TRUE)
clu <- c(0.3)
#methods <- c("cca", "harmony", "rpca")
methods <- c( "harmony")
#split <- "Sample"
proj.name <- "AFL_NEW_All"
specie<- "hsa" # homo sapiens "hsa" or mus musculus "mmu"
integrated_file_name <- "AFL_NEW_All.integrated_od_markers_annotated_09.05.2025.rds" 
################# SETTING VARIABLES ##############################################
#added the last row to the color vector
custom_colors <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51',
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62',
  '#e8e0d8','#fcf8f2','#eae4b7','#79958a','#345146'
)
#######################    FUNCTIONS      ######################################
PlotBarplots <- function(data, x, group, plottype = c("fill","stack","dodge"), colvec){
  
  dat <- aggregate(list("Frequency" = data[[group]]), data[c(group,x)], FUN = length)
  
  p2.data <- ggplot(dat, aes(fill=Sample, y=Frequency, x=integrated_snn_res.0.3)) +
    geom_bar(position=plottype, stat="identity") +
    scale_fill_manual("legend", values = colvec)
  p2.data
}
save_plot <- function(path=output.path, plot=NULL, plot_name=NULL){
  plot_name <- paste0(paste(plot_name,today,sep = "_"), ".pdf")
  ggsave(filename= plot_name , plot=plot, path=path, width = 14, height = 8, dpi = 300) 
}
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)
}
#######################    FUNCTIONS      ######################################
#################    READ SEURAT LIST     ######################################
#Read seurat object list all samples
integrated <- readRDS(paste(srcRaw, integrated_file_name, sep="/"))
#################    READ SEURAT LIST     ######################################

################# CHECK CELL DISTRIBUTION TO SELECT USED TEST  #################
# Extract the raw counts matrix from your Seurat object
counts_matrix <- GetAssayData(object = integrated, slot = "counts")
# Flatten the counts matrix into a vector
counts_vector <- as.vector(counts_matrix)
# Plot a histogram of raw counts
hist(counts_vector, breaks = 100, main = "Histogram of Raw Counts", 
     xlab = "Counts", ylab = "Frequency", col = "skyblue", border = "white")
# Extract the normalized data (log-normalized counts)
normalized_matrix <- GetAssayData(object = integrated, slot = "data")  # "data" slot contains log-normalized values
normalized_vector <- as.vector(normalized_matrix)
# Plot a histogram of log-normalized counts
hist(normalized_vector, breaks = 100, main = "Histogram of Log-Normalized Counts", 
     xlab = "Log-Normalized Counts", ylab = "Frequency", col = "lightgreen", border = "white")
# Calculate the proportion of zero counts
zero_fraction <- sum(counts_matrix == 0) / length(counts_matrix)
cat("Proportion of zero counts:", zero_fraction, "\n")
#Proportion of zero counts: 0.8306797 
#A proportion of zero counts of ~83% indicates that your data is highly zero-inflated, 
#which is typical for single-cell RNA-seq data. 
#MAST - Best Choice for Zero-inflated Data
#Reasons to use MAST:
#  Zero-inflation handling:
#MAST uses a hurdle model that separately accounts for the presence/absence of expression (binary, zero vs non-zero) and the actual expression level (continuous). This is ideal for your data with 86% zeros.
#Latent variable regression:
#MAST allows you to include confounding factors (latent.vars), such as sex, which you want to regress out.
#Widely used for scRNA-seq:
#MAST is specifically designed for single-cell data and is one of the most robust and commonly used methods in the field.
################# CHECK CELL DISTRIBUTION TO SELECT USED TEST  #################



########################      SINGLE CELL DGE ANALYSIS ########################      
Idents(integrated) <- "integrated.clusters.harmony.Batch"
FindMarkersMAST <- list()
#Comparisons to make, alwasy will take the first as reference
comparisons <- integrated$Condition %>% unique(.)
ref <- comparisons[1]
comparisons <- comparisons[-1]
for (comp in comparisons){
  FindMarkersMAST[[paste0(ref,"_vs_",comp)]] <- FindMarkers(object = integrated, 
                                                            group.by= "Condition",
                                                            ident.1 = ref,#LogFC Possitive
                                                            ident.2 = comp, #Reference
                                                            test.use = "MAST",
                                                            assay="RNA",
                                                            #latent.vars = c("Batch")
  )
  FindMarkersMAST[[paste0(ref,"_vs_",comp)]]$genes <- rownames(FindMarkersMAST[[paste0(ref,"_vs_",comp)]])
}
#Save the find markers 
print("Saving FindMarkers into an excel (.xlsx) file...")
#list_of_datasets <- list("Condition_KOvsWT" = FindMarkersMAST)
writexl::write_xlsx(FindMarkersMAST,
                    paste(output.path,"FindMarkersMAST_Condition_comparison.xlsx",sep = "/"))
########################     VISUALIZATION        ##############################
library(EnhancedVolcano)
volcano_plot_pvalue <- list()
volcano_plot_Adjpvalue <- list()
for (comp in comparisons){
  table <- FindMarkersMAST[[paste0(ref,"_vs_",comp)]]
  volcano_plot_pvalue[[paste0(ref,"_vs_",comp)]]<-EnhancedVolcano(table,
                                                                  lab = table$genes,
                                                                  x = 'avg_log2FC',
                                                                  y = 'p_val',
                                                                  title = paste0(ref," vs ",comp),
                                                                  subtitle= paste0("p-val: 0.05"),
                                                                  pCutoff = 0.05,
                                                                  #ylim=c(0,300)
  )
  volcano_plot_Adjpvalue[[paste0(ref,"_vs_",comp)]]<-EnhancedVolcano(table,
                                                                     lab = table$genes,
                                                                     x = 'avg_log2FC',
                                                                     y = 'p_val_adj',
                                                                     title = paste0(ref," vs ",comp),
                                                                     subtitle= paste0("Adjusted p-val: 0.05"),
                                                                     pCutoff = 0.05,
                                                                     #ylim=c(0,300)
  )
}
layout <- "1122"
plot_name <- paste(output.path,"scRNAseq_Conditions_Volcano.pdf", sep = "/")
pdf(file= plot_name,onefile = TRUE, width = 18, height = 8)
for (comp in comparisons){
  caca<- volcano_plot_pvalue[[paste0(ref,"_vs_",comp)]] + volcano_plot_Adjpvalue[[paste0(ref,"_vs_",comp)]] +plot_layout(design = layout)+
    plot_annotation(caption= paste0("Date of render: ",today, "\n"))
  print(caca)
}
dev.off()
########################     VISUALIZATION        ##############################
########################      SINGLE CELL DGE ANALYSIS WITHIN EACH CLUSTER ########################      


########################     VISUALIZATION PCA CORRECTION       ################
sample_PCA_bef <- DimPlot(integrated, reduction ="pca", group.by = "Sample" )+labs(title = "PCA No Correction")
sample_PCA_corr <- DimPlot(integrated, reduction ="integrated.harmony.Batch", group.by = "Sample" )+labs(title = "PCA Corrected")
batch_PCA_bef <- DimPlot(integrated, reduction ="pca", group.by = "Batch" )+labs(title = "PCA No Correction")
batch_PCA_corr <- DimPlot(integrated, reduction ="integrated.harmony.Batch", group.by = "Batch" )+labs(title = "PCA Corrected")
p <-  sample_PCA_bef | sample_PCA_corr
save_plot(path=output.path, plot=p, plot_name="PCA_SampleComparison.pdf")
########################     VISUALIZATION PCA CORRECTION       ################