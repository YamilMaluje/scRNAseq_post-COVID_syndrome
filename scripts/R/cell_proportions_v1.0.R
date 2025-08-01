
#TUTORIAL 
#https://phipsonlab.github.io/propeller-paper-analysis/RealDataAnalysis.html#COVID_data
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
#Propeller packages
library(speckle)
library(limma)
library(edgeR)
library(pheatmap)
library(gt)
#scproportiontest package
library("scProportionTest")
################# SETTING VARIABLES ##############################################
orig <- getwd()
today <- Sys.Date() %>% format(., format="%d.%m.%Y")
output.path <- paste(orig,"05_DownstreamAnalysis","Diff_CellType_Population",sep = "/")
output.path.analysis <- paste(output.path, today,sep = "/")
#Tool of analysis 
#output.path.analysis.tool <- paste(output.path.analysis, "propeller_annot",sep = "/")
output.path.analysis.tool <- paste(output.path.analysis, "sc_proportion_annot",sep = "/")
dir.create(output.path.analysis.tool,recursive = TRUE)
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
save_plot <- function(path=output.path, plot=NULL, plot_name=NULL){
  plot_name <- paste0(paste(plot_name,today,sep = "_"), ".pdf")
  ggsave(filename= plot_name , plot=plot, path=path, width = 14, height = 7, dpi = 200) 
}
#######################    FUNCTIONS      ######################################
#################    READ SEURAT LIST     ######################################
obj <- readRDS("/Volumes/MarkV/Projects/UKSH scRNAseq/Server Files/scRNAseq Analysis/COVID/COVID-AFL_NEW_ALL/05_DownstreamAnalysis/SelectedGenes/09.05.2025/AFL_NEW_All.integrated_od_markers_annotated_09.05.2025.rds")
obj$Sample_fac <- factor(obj$Sample)
obj$Condition_fac <- factor(obj$Condition)
#################    READ SEURAT LIST     ######################################

iden <- "CanonicalMarkers_Annotation"
########################    MAIN CODE     ######################################
# PROPELLER   #############################
Idents(obj) <- iden
output.logit <- propeller(x=obj, clusters=factor(obj@meta.data[[iden]]), sample=obj$Sample_fac, group = obj$Condition_fac, transform="asin")

#Following this tutorial
#https://phipsonlab.github.io/propeller-paper-analysis/pbmcJP.html
#Explore cell type proportions among the Condition 
props <- getTransformedProps(clusters = obj@meta.data[[iden]], 
                             sample = obj$Condition)
p1 <- plotCellTypeProps(clusters = obj@meta.data[[iden]], sample = obj$Condition) + 
  theme(axis.text.x = element_text(angle = 45))+ 
  ggtitle("Cell Cluster proportions") + 
  theme(plot.title = element_text(size = 18, hjust = 0))
p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = 0))
save_plot(path = output.path.analysis.tool, plot_name = "Propeller_CellCluster_proportion", plot = p1)
#Exploring heterogeneity in cell type proportions between Condition
counts <- table(obj@meta.data[[iden]], obj$Condition)
baselineN <- rowSums(counts)
N <- sum(baselineN)
baselineprops <- baselineN/N
props <- getTransformedProps(clusters = obj@meta.data[[iden]], 
                             sample = obj$Condition)
cols <- ggplotColors(nrow(props$Proportions))
m <- match(rownames(props$Proportions),levels(factor(obj@meta.data[[iden]])))
pdf(paste(output.path.analysis.tool, "Propeller_CellCluster_proportion_jitter.pdf", sep="/"))
par(mfrow=c(1,1))
par(mar=c(7,5,2,2))
plot(jitter(props$Proportions[,1]), col = cols[m], pch=16, ylim=c(0,max(props$Proportions)),
     xaxt="n", xlab="", ylab="Cell cluster proportion", cex.lab=1.5, cex.axis=1.5)
for(i in 2:ncol(props$Proportions)){
  points(jitter(1:nrow(props$Proportions)),props$Proportions[,i], col = cols[m],
         pch=16)
}
axis(side=1, at=1:nrow(props$Proportions), las=2, 
     labels=rownames(props$Proportions))
title("Cell cluster proportions estimates for 4 Conditions")
dev.off()
# data is overdispersed compared to what would be expected under a Binomial or Poisson distribution.
pdf(paste(output.path.analysis.tool, "Propeller_mean-variance_relationship_counts.pdf", sep="/"))
plotCellTypeMeanVar(counts)
dev.off()
pdf(paste(output.path.analysis.tool, "Propeller_mean-variance_relationship_counts_proportion.pdf", sep="/"))
plotCellTypePropsMeanVar(counts)
dev.off()



df_list<- list()
df_list$Counts <-as.data.frame.matrix(props[["Counts"]]) 
df_list$Counts["Cluster"] <- rownames(df_list$Counts)
df_list$TransformedProps <-as.data.frame.matrix(props[["TransformedProps"]])
df_list$TransformedProps["Cluster"] <- rownames(df_list$TransformedProps)
df_list$Proportions <-as.data.frame.matrix(props[["Proportions"]])
df_list$Proportions["Cluster"] <- rownames(df_list$Proportions)
slot.name <- "Propeller_proportion"
writexl::write_xlsx(df_list,path = paste(output.path.analysis.tool,
                                         paste0(slot.name,today,".xlsx"),
                                         sep = "/"),
)






#SC PROPORTION  ###########################
#https://github.com/rpolicastro/scProportionTest
#Explanation
#This R library facilitates the analysis of the difference between the proprotion of cells in clusters between two scRNA-seq samples. 
#A permutation test is used to calculate a p-value for each cluster, and a confidence interval for the magnitude difference is returned via bootstrapping. 
table_list <- list()
prop_test <- sc_utils(obj)
table(prop_test@meta_data[["Condition"]])
#PBS      TGFb      TNFa TNFa_TGFb 
#7942      4059      8707      2247 
#Comparisons to make, alwasy will take the first as reference
comparisons <- obj$Condition %>% unique(.)
ref <- comparisons[1]
comparisons <- comparisons[-1]
for(comp in comparisons){
  prop_test <- permutation_test(prop_test, 
                                cluster_identity = iden,
                                sample_1 = ref, 
                                sample_2 = comp,
                                sample_identity = "Condition",
                                n_permutations = 1000)
  
  plot  <- permutation_plot(prop_test,log2FD_threshold = log2(1.415)) + 
    ggtitle(label = "scPropotionTest",subtitle = paste(comp,ref,sep = " v/s ")) +
    labs(caption= paste0("Date of render: ",today))
  save_plot(path = output.path.analysis.tool, plot_name = paste0("scProportion_",comp,"_vs_",ref), plot = plot)
  
  table_list[[paste(comp,ref,sep = "_vs_")]] <- prop_test@results[["permutation"]]
  
}

slot.name <- "scProportion_comparison_"
writexl::write_xlsx(table_list,path = paste(output.path.analysis.tool,
                                            paste0(slot.name,today,".xlsx"),
                                            sep = "/"))
########################    MAIN CODE     ######################################
