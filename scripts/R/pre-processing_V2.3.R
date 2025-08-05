################# SCRIPT INFORMATION  ##############################################
#This script will read the count matrices and perform BarcodeRanks or DropletUtils::emptyDrops depending
#on the sequencing technology (10x, BD, Singleron, ...) to remove the empty droplets
#Then will identify possible doublets with "Main_Doublet_Decon()" and tag each cell as possible doublet or not
#Then will remove the doublets/multiplet from the seurat object
#Then perform the normal seurat workflow:
#PercentageFeatureSet()
#NormalizeData
#FindVariableFeatures
#ScaleData
#RunPCA
#RunTSNE
#RunUMAP
#FindNeighbors
#FindClusters
#Then will subset the seurat object depending on the selected mitochondrial percentage, Number of counts and number of features from an Excel file were the sheet must be called "QC"
#Along each process the number of cells is stored as well as different box plots, violin plots, scatter plots are created for debuging 
#Also a pdf and xlsx file is created with the the number of cells along each step and also a overview of the entire dataset percentage is calculated to quantify
#the percentage of lost cells compared with the starting point 
#NOTE: BEFORE RUNING THIS YOU NEED TO SELECT THE SUBSET METHOD TO USE MANUALLY OR AUTOMATIC
################# SCRIPT INFORMATION  ##############################################
library(Seurat)
library(dplyr)
library(Matrix)
library(DropletUtils)
library(tidyverse)
library(ggpointdensity)
library(scico)
library(scales)
library(BUSpaRse)
#library(scater)
library(ggplot2)
library(patchwork)
library(scDblFinder)
library(gridExtra)
library(openxlsx)
library(ggpubr)  # For fitting the curve
################# SETTING VARIABLES ##############################################
orig <- getwd()
srcRaw <- "02_align"
# Define the patterns for folder names
#pattern1 <- "(C25010700[0-9]{1})_scR"
pattern1 <- "(C250[0-9]{6})_scR"

pattern2 <- "cells_x_genes"
# List all files and directories in the specified path
all_paths <- list.dirs(orig, recursive = TRUE, full.names = TRUE)
# Filter paths based on the specified conditions
filtered_paths <- grep(paste0("/", pattern1, "/.*", pattern2), all_paths, value = TRUE)
# Make the paths unique
mypath <- unique(filtered_paths)
today <- Sys.Date() %>% format(., format="%d.%m.%Y")
#Create output folders to save all the files
output.path <- paste(orig,"03_QC",today ,sep = "/")
output.path.plots <- paste(output.path,"Plots",sep = "/")
dir.create(output.path.plots,recursive = TRUE)
proj.name <- "COVID-AFL_NEW_All"
specie<- "hsa" # homo sapiens "hsa" or mus musculus "mmu"
#gns <- read.delim(file = paste(orig,"tr2g_hgnc.MM.txt",sep = "/"), header = TRUE, col.names = c("transcript","gene"))
#gns <- read.delim(file = paste(orig,"t2g.txt",sep = "/"), header = F, col.names = c("transcript","gene","gene_name"))
seq_tech="Singleron" #which sequencing technology was used "10x","BD","Singleron"....
#Metadata file where in the first sheet is the metadata and second sheet "QC" are the QC parameters to set per individual sample
metadata_file <- "metadata_covid_TNF_NEW_All.xlsx"
################# SETTING VARIABLES ##############################################
################# FUNCTIONS ##############################################
knee_plot <- function(bc_ranks) {
  knee_plt <- tibble(rank = map(bc_ranks, ~ .x[["rank"]]), 
                     total = map(bc_ranks, ~ .x[["total"]]),
                     dataset = names(bc_ranks)) %>% 
    unnest(cols = c(rank, total)) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = map_dbl(bc_ranks, ~ metadata(.x)[["inflection"]]),
                  rank_cutoff = map_dbl(bc_ranks, 
                                        ~ max(.x$rank[.x$total >
                                                        metadata(.x)[["inflection"]]])),
                  dataset = names(bc_ranks))
  p <- ggplot(knee_plt, aes(rank, total, color = dataset)) +
    geom_line() +
    geom_hline(aes(yintercept = inflection, color = dataset), 
               data = annot, linetype = 2) +
    geom_vline(aes(xintercept = rank_cutoff, color = dataset),
               data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Rank", y = "Total UMIs")
  return(p)
}
seu_plots <- function(seu.obj.list, output.path =".", processing.stage=""){
  library(patchwork)
  genes.boxplot <- data.frame() 
  for(i in 1:length(seu.obj.list)){
    aux <- seu.obj.list[[i]][[c("Sample","percent.mt", "nFeature_RNA","nCount_RNA")]] %>% as.data.frame()
    genes.boxplot <- rbind(genes.boxplot,aux)
  }
  
  # COUNT NUMBER
  p1 <- ggplot(genes.boxplot, aes(x=Sample, y=nCount_RNA, fill=Sample)) + 
    geom_boxplot() +
    ggtitle("Count Number", subtitle = "All samples") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  # MT percentage
  p2 <- ggplot(genes.boxplot, aes(x=Sample, y=percent.mt, fill=Sample)) + 
    geom_boxplot() +
    ggtitle("Mitochondria percentage", subtitle = "All samples") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  # Number Features
  p3 <- ggplot(genes.boxplot, aes(x=Sample, y=nFeature_RNA, fill=Sample)) + 
    geom_boxplot() +
    ggtitle("Feature Number", subtitle = "All samples") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  name <- paste0(processing.stage, "_BoxPlot_Number_Cell.MT.Features_allSamples.pdf")
  pdf(paste(output.path, name , sep = "/" ))
  print(list(p1,p3,p2))
  dev.off()
  
  #Violin plot  OF ALL THE SAMPLES
  #This to check the mt threshold
  v_plot <-list()
  for(i in 1:length(seu.obj.list)){
    Idents(seu.obj.list[[i]]) <- "orig.ident"
    p <- VlnPlot(object= seu.obj.list[[i]], features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "orig.ident") +
      plot_annotation(title = paste0("Sample: ",names(seu.obj.list)[i]))
    
    #Scatter plot
    plot1 <- FeatureScatter(seu.obj.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    plot2 <- FeatureScatter(seu.obj.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    plot <- CombinePlots(plots = list(plot1, plot2))
    
    #equals to Feature Scatter function but with MT color code
    plot3 <- ggplot(seu.obj.list[[i]]@meta.data, aes(nCount_RNA, nFeature_RNA)) +   #### Repeat for the other 3 samples (10x$FVB, BD$Bl6, 10x$FVB)  
      geom_point(alpha = 0.7, size = 0.5, aes(colour = percent.mt)) +
      scale_colour_gradient(low = "blue", high = "yellow") +
      labs(x = "Total UMI counts per cell", y = "Number of genes detected") + 
      ggtitle(label = paste0("Genes to Transcripts Scatterplot ",names(seu.obj.list[[i]])), subtitle = paste0("Sample: ",names(seu.obj.list)[i]))
    
    
    v_plot[[names(seu.obj.list)[i] ]] <- list("ViolinPlot"=p,
                                              "ScatterPlot"=plot,
                                              "ScatterPlot2"=plot3)
    
    
  }
  name <- paste0(processing.stage, "_ViolinPlot_Number_Cell.MT.Features_allSamples.pdf")
  pdf(paste(output.path, name , sep = "/" ))
  print(v_plot)
  dev.off()
  rm(v_plot)
}
empty_drops_removal <- function (count_mat=NULL,sample="Sample"){
  #count_mat= matrix containing UMI counts for all barcodes
  npts= 10000
  set.seed(100)
  e.out <- DropletUtils::emptyDrops(count_mat,niters = npts)
  
  is.cell <- e.out$FDR <= 0.01
  sumcells<- sum(is.cell, na.rm=TRUE)
  cat("\nEmpty Droplets sample: ",sample)
  cat("\nNumber of possible cells:", sumcells)
  
  # The Limited field indicates whether a lower p-value could be obtained by increasing the number of permutations.
  #If there are any entries with FDR above the desired threshold and Limited==TRUE, it indicates that npts should be increased in the emptyDrops call.
  df <-as.data.frame(table(Limited=e.out$Limited, Significant=is.cell))
  
  
  value <- df[df$Significant == FALSE & df$Limited == TRUE, ]$Freq 
  while(value > 100){
    #we have to change the number of number of permutation ins empty dropsnpts
    npts= as.integer(npts*1.10)
    cat("\n ingreasing the number of permutations for the Monte Carlo p-value calculation to ", npts)
    e.out <- DropletUtils::emptyDrops(count_mat,niters = npts)
    
    is.cell <- e.out$FDR <= 0.01
    sumcells<- sum(is.cell, na.rm=TRUE)
    
    cat("\nEmpty Droplets sample: ",sample)
    cat("\nNumber of possible cells:", sumcells)
    df <-as.data.frame(table(Limited=e.out$Limited, Significant=is.cell))
    
    value <- df[df$Significant == FALSE & df$Limited == TRUE, ]$Freq 
    
  }
  
  #Droplets detected as cells should show up with large negative log-probabilities or very large total counts
  df <- data.frame(Total = e.out$Total,
                   LogProb = -e.out$LogProb,
                   is.cell = is.cell)
  p <- ggplot(data= df, aes(x = Total, y = LogProb, color = is.cell)) +
    geom_point() +
    scale_color_manual(values = c("black", "red")) +
    labs(x = "Total UMI count",
         y = "-Log Probability",
         title = paste0("Sample: ", sample),
         caption = paste0("# of Cells: ", sumcells))
  
  
  #remove empty droplets
  e.out.iscell <- subset(e.out, subset= FDR <=0.01) %>% rownames(.)
  count_mat_iscell <- count_mat[, colnames(count_mat) %in% e.out.iscell] 
  l <- list("plot"=p,
            "count_mat_iscell"=count_mat_iscell )
  
  return(l)
}
seu_plots_withThresholdLines <- function(seu.obj,qc_parameter_threshold){
  library(patchwork)
  
  #Violin plot  
  Idents(seu.obj) <- "orig.ident"
  p1 <- VlnPlot(object= seu.obj, features = c("nCount_RNA"), group.by = "orig.ident") +
    plot_annotation(title = paste0("Sample: ",seu.obj$Sample[1])) +
    # Add vertical and horizontal lines at the inflection point
    geom_hline(yintercept = qc_parameter_threshold$nCount_RNA_Low, linetype = "dashed", color = "red") +
    geom_hline(yintercept = qc_parameter_threshold$nCount_RNA_High, linetype = "dashed", color = "red") +
    NoLegend()
  
  p2 <- VlnPlot(object= seu.obj, features = c("nFeature_RNA"), group.by = "orig.ident") +
    plot_annotation(title = " ") +
    # Add vertical and horizontal lines at the inflection point
    geom_hline(yintercept = qc_parameter_threshold$nFeature_RNA_Low, linetype = "dashed", color = "red") +
    geom_hline(yintercept = qc_parameter_threshold$nFeature_RNA_High, linetype = "dashed", color = "red") +
    NoLegend()
  
  p3 <- VlnPlot(object= seu.obj, features = c("percent.mt"), group.by = "orig.ident") +
    plot_annotation(title = " ") +
    # Add vertical and horizontal lines at the inflection point
    geom_hline(yintercept = qc_parameter_threshold$MT_High, linetype = "dashed", color = "red") +
    NoLegend()
  p <- CombinePlots(plots = list(p1, p2, p3),ncol = 3)
  
  #Scatter plot
  plot1 <- FeatureScatter(seu.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    NoLegend()
  plot2 <- FeatureScatter(seu.obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    NoLegend()
  plot <- CombinePlots(plots = list(plot1, plot2))
  
  #equals to Feature Scatter function but with MT color code
  plot3 <- ggplot(seu.obj@meta.data, aes(nCount_RNA, nFeature_RNA)) +   #### Repeat for the other 3 samples (10x$FVB, BD$Bl6, 10x$FVB)  
    geom_point(alpha = 0.7, size = 0.5, aes(colour = percent.mt)) +
    scale_colour_gradient(low = "blue", high = "yellow") +
    labs(x = "Total UMI counts per cell", y = "Number of genes detected") + 
    ggtitle(label = "Genes to Transcripts Scatterplot ", 
            subtitle = paste0("Sample: ",seu.obj$Sample[1]))+
    geom_vline(xintercept = qc_parameter_threshold$nCount_RNA_Low, linetype = "dashed", color = "red") +
    geom_vline(xintercept = qc_parameter_threshold$nCount_RNA_High, linetype = "dashed", color = "red") +
    geom_hline(yintercept = qc_parameter_threshold$nFeature_RNA_Low, linetype = "dashed", color = "red") +
    geom_hline(yintercept = qc_parameter_threshold$nFeature_RNA_High, linetype = "dashed", color = "red") 
  
  v_plot <- list("ViolinPlot"=p,
                 "ScatterPlot"=plot,
                 "ScatterPlot2"=plot3)
  return (v_plot)
}
seu_plots_withThresholdLines_method1 <- function(seu.obj,qc_parameter_threshold){
  library(patchwork)
  
  #Violin plot  
  Idents(seu.obj) <- "orig.ident"
  p1 <- VlnPlot(object= seu.obj, features = c("nCount_RNA"), group.by = "orig.ident") +
    plot_annotation(title = paste0("Sample: ",seu.obj$Sample[1])) +
    # Add vertical and horizontal lines at the inflection point
    geom_hline(yintercept = qc_parameter_threshold$nCount_RNA_High, linetype = "dashed", color = "red") +
    NoLegend()
  
  p2 <- VlnPlot(object= seu.obj, features = c("nFeature_RNA"), group.by = "orig.ident") +
    plot_annotation(title = " ") +
    # Add vertical and horizontal lines at the inflection point
    geom_hline(yintercept = qc_parameter_threshold$nFeature_RNA_High, linetype = "dashed", color = "red") +
    NoLegend()
  
  p3 <- VlnPlot(object= seu.obj, features = c("percent.mt"), group.by = "orig.ident") +
    plot_annotation(title = " ") +
    # Add vertical and horizontal lines at the inflection point
    NoLegend()
  p <- CombinePlots(plots = list(p1, p2, p3),ncol = 3)
  
  #Scatter plot
  plot1 <- FeatureScatter(seu.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    NoLegend()
  plot2 <- FeatureScatter(seu.obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    NoLegend()
  plot <- CombinePlots(plots = list(plot1, plot2))
  
  #equals to Feature Scatter function but with MT color code
  plot3 <- ggplot(seu.obj@meta.data, aes(nCount_RNA, nFeature_RNA)) +   #### Repeat for the other 3 samples (10x$FVB, BD$Bl6, 10x$FVB)  
    geom_point(alpha = 0.7, size = 0.5, aes(colour = percent.mt)) +
    scale_colour_gradient(low = "blue", high = "yellow") +
    labs(x = "Total UMI counts per cell", y = "Number of genes detected") + 
    ggtitle(label = "Genes to Transcripts Scatterplot ", 
            subtitle = paste0("Sample: ",seu.obj$Sample[1]))+
    geom_vline(xintercept = qc_parameter_threshold$nCount_RNA_High, linetype = "dashed", color = "red") +
    geom_hline(yintercept = qc_parameter_threshold$nFeature_RNA_High, linetype = "dashed", color = "red") 
  
  v_plot <- list("ViolinPlot"=p,
                 "ScatterPlot"=plot,
                 "ScatterPlot2"=plot3)
  return (v_plot)
}
analyze_scatterplot <- function(seurat_object,probs=0.1,model_name="polyreg") {
  library(stats)
  
  # Create a scatter plot
  p <- ggplot(seurat_object@meta.data, aes(nCount_RNA, nFeature_RNA)) +
    geom_point(alpha = 0.7, size = 0.5, aes(colour = percent.mt)) +
    scale_colour_gradient(low = "blue", high = "yellow") +
    labs(x = "Total UMI counts per cell", y = "Number of genes detected") +
    ggtitle(label = paste0("Genes to Transcripts Scatterplot ", names(seurat_object)), 
            subtitle = paste0("Sample: ", seurat_object$Sample[1]))
  # Calculate a non-linear fit
  model <-switch (model_name,
                  #Determine the nonlinear (weighted) least-squares estimates of the parameters of a nonlinear model.
                  "nls" = nls(nFeature_RNA ~ a * nCount_RNA^b, data = seurat_object@meta.data, start = list(a = 1, b = 1)),
                  #Polynomial regression fits a polynomial function to the data. We'll fit a quadratic polynomial in this example. 
                  "polyreg"= lm(nFeature_RNA ~ poly(nCount_RNA, 3), data = seurat_object@meta.data),
                  #linear model
                  "lm" = lm(nFeature_RNA ~ poly(nCount_RNA, 1), data = seurat_object@meta.data),
                  
  )
  # Predicted values
  seurat_object@meta.data$predicted <- predict(model)
  
  # Calculate residuals
  seurat_object@meta.data$residuals <- residuals(model)
  
  # Calculate threshold
  threshold <- quantile(seurat_object@meta.data$residuals, probs = probs)
  # Create a new column for viability
  seurat_object@meta.data$viable <- ifelse(abs(seurat_object@meta.data$residuals) <= abs(threshold), "Viable", "Non-Viable")
  table(seurat_object@meta.data$viable)
  
  # Filter cells above threshold
  viable_cells <- seurat_object@meta.data %>%
    filter(viable == "Viable")
  
  # Find the point for the threshold
  # Find the cell with the biggest nCount_RNA and nFeature_RNA
  threshold_point <- viable_cells %>%
    slice_max(order_by = nCount_RNA, n = 1) %>%
    slice_max(order_by = nFeature_RNA, n = 1)%>%
    select(nCount_RNA, nFeature_RNA)
  threshold_point$Sample <-  seurat_object@meta.data$Sample[1]
  # Add the fitted curve
  curve_df <- data.frame(nCount_RNA = seq(min(seurat_object@meta.data$nCount_RNA), max(seurat_object@meta.data$nCount_RNA), length.out = 100))
  curve_df$nFeature_RNA <- predict(model, newdata = curve_df)
  
  #This plot is only to check if it is correct
  #p <- ggplot(seurat_object@meta.data, aes(nCount_RNA, nFeature_RNA)) +
  #  geom_point(aes(colour = viable), alpha = 0.7, size = 0.5) +
  #  scale_colour_manual(values = c("Viable" = "blue", "Non-Viable" = "grey"), guide = FALSE) +
  #  labs(x = "Total UMI counts per cell", y = "Number of genes detected") +
  #  ggtitle(label = paste0("Genes to Transcripts Scatterplot ", names(seurat_object)), 
  #          subtitle = paste0("Sample: ", seurat_object$Sample[1]))
  p <- p +
    geom_line(data = curve_df, aes(nCount_RNA, nFeature_RNA), color = "black", size = 1) +  # Fitted curve
    geom_vline(xintercept = threshold_point$nCount_RNA, linetype = "dashed") +
    geom_hline(yintercept = threshold_point$nFeature_RNA, linetype = "dashed")
  
  out <- list("plot"=p,
              "threshold"= threshold_point,
              "slected_cells"= viable_cells)
  return(out)
}
################# FUNCTIONS ##############################################
#loading indepth tr2g file as lookup table
#ah <- AnnotationHub(cache = paste0("/Users/fabianott/work/SingleCell/Rstuff/AnnotationHub/"))
#query(ah, pattern = c("Ensembl", "104", "Homo Sapiens", "EnsDb"))
#edb <- ah[["AH75011"]]
#gns <- tr2g_EnsDb(ensdb = edb,Genome = BSgenome.Hsapiens.UCSC.hg38,use_gene_version = FALSE)[,c("transcript","gene", "gene_name")] %>% distinct()
#appendSample <- F
#if (appendSample == T){
#  #MatchSamples <- c("C33/cells_x_genes","C34/cells_x_genes","C35/cells_x_genes","C36/cells_x_genes","C37/cells_x_genes","C38/cells_x_genes")
#  #MatchSamples <- c("C33/cells_x_genes","C34/cells_x_genes")
#  MatchSamples <- c("S7_4/Solo.out/Gene/raw","S7_12/Solo.out/Gene/raw","S7_20/Solo.out/Gene/raw","S7_28/Solo.out/Gene/raw")
#  
#  toMatch <- MatchSamples %>% paste(. , collapse = "|")
#  mypath <- mypath[grep(pattern = toMatch, x = mypath)]
#  #read existing seurat object list
#  obj.list <- readRDS(paste(orig,"noDub.covList.20220208.rds",sep = "/"))
#} else {
#  seu.list <- list()
#}
###################   READ COUNT MATRIX   ######################
seu.list <- list()
res_mat.raw <- list()
res_mat <- list()
plot.BCRanks <- list()
#i=mypath[[1]]
for (i in mypath){
  #obtain the folder name that is equivalent to the sample id
  nm <- str_extract(i, pattern = pattern1)
  message(sprintf("Processing sample= %s",nm))
  
  #which sequencing technology was used
  switch(seq_tech,
         "10x" = {
           
           res_mat.raw[[nm]] <- Read10X(i,
                                        gene.column = 2,
                                        cell.column = 1,  
                                        unique.features = TRUE,
                                        strip.suffix = FALSE)
           
           
           temp <- empty_drops_removal(count_mat=res_mat.raw[[nm]], sample=nm )
           
           plot.BCRanks[[nm]] <- temp$plot
           res_mat[[nm]] <- temp$count_mat_iscell
           
         },
         
         "BD"= { 
           i.temp <- paste(i,"cells_x_genes",sep = "/")
           # Read in sparse matrix transcripts
           res_mat.raw[[nm]] <- read_count_output(dir = i.temp,
                                                  name = "",
                                                  tcc = FALSE)
           
           
           # Inflection point hard coded
           totalCounts <- Matrix::colSums(res_mat.raw[[nm]] )
           bc_rank <- barcodeRanks(res_mat.raw[[nm]] ) 
           # Kneeplot (ggplot style)
           plot.BCRanks[[nm]] <- qplot(bc_rank$total, bc_rank$rank, geom= "line") +
             geom_vline(xintercept = metadata(bc_rank)$knee, color="blue", linetype=2) +
             geom_vline(xintercept = metadata(bc_rank)$inflection, color = "green", linetype=2) +
             annotate("text", y=1000, x=1.5 *c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection),
                      label=c("knee", "inflection"), color = c("blue","green")) +
             scale_x_log10() + scale_y_log10() +
             labs(y="Barcode rank", x = "Total UMI count",title = paste0("Sample: ", nm))
           
           #########   MOD SECTION FOR MANUAL SELECT INFLECTION POINT
           #infl.point <- 2000
           #plot.BCRanks[[nm]] <- qplot(bc_rank$total, bc_rank$rank, geom= "line") +
           #  geom_vline(xintercept = metadata(bc_rank)$knee, color="blue", linetype=2) +
           #  geom_vline(xintercept = infl.point, color = "green", linetype=2) +
           #  annotate("text", y=1000, x=1.5 *c(metadata(bc_rank)$knee, infl.point),
           #           label=c("knee", "inflection"), color = c("blue","green")) +
           #  scale_x_log10() + scale_y_log10() +
           #  labs(y="Barcode rank", x = "Total UMI count",title = paste0("Sample: ", nm))
           
           #res_mat[[nm]] <- res_mat.raw[[nm]][,totalCounts > infl.point]
           #######   MOD SECTION FOR MANUAL SELECT INFLECTION POINT
           
           res_mat[[nm]] <- res_mat.raw[[nm]][,totalCounts > bc_rank@metadata$inflection]
           
         },
         
         
         "Singleron"={
           gns <- read.delim(file = paste(orig,"tr2g_human.tsv",sep = "/"), header = TRUE, col.names = c("transcript","gene","gene_name"))
           
           # Read in sparse matrix transcripts
           res_mat.raw[[nm]] <- read_count_output(dir = i,
                                                  name = "cells_x_genes",
                                                  tcc = FALSE)
           
           # Do some HGNC reformating stuff
           chck1 <- dim(res_mat.raw[[nm]])
           idx <- rownames(res_mat.raw[[nm]])
           idx <- as.data.frame(x = idx, stringsAsFactors=FALSE)
           names(idx) <- "gene"
           SymID <- plyr::join(x = idx, y = gns, by = "gene")
           SymID[SymID$gene_name %>% is.na(.), 3 ] <- SymID[SymID$gene_name %>% is.na(.), 1 ]
           SymID <- SymID[!duplicated(x = SymID$gene),]
           SymID[SymID$gene_name %>% '=='(.,""), 3] <- SymID[SymID$gene_name %>% '=='(.,""), 1]
           SymID[SymID$gene_name %>% duplicated(.),3] <- paste0(SymID[SymID$gene_name %>% duplicated(.),3],"-1")
           rownames(res_mat.raw[[nm]]) <- SymID$gene_name
           chck2 <- dim(res_mat.raw[[nm]])
           if(chck1[1] == chck2[1]){
             message(sprintf("Replacing Ensembl with HGNC went well!"))
           }
           
           # Inflection point hard coded
           totalCounts <- Matrix::colSums(res_mat.raw[[nm]] )
           bc_rank <- barcodeRanks(res_mat.raw[[nm]] ) 
           # Kneeplot (ggplot style)
           plot.BCRanks[[nm]] <- qplot(bc_rank$total, bc_rank$rank, geom= "line") +
             geom_vline(xintercept = metadata(bc_rank)$knee, color="blue", linetype=2) +
             geom_vline(xintercept = metadata(bc_rank)$inflection, color = "green", linetype=2) +
             annotate("text", y=1000, x=1.5 *c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection),
                      label=c("knee", "inflection"), color = c("blue","green")) +
             scale_x_log10() + scale_y_log10() +
             labs(y="Barcode rank", x = "Total UMI count",title = paste0("Sample: ", nm))
           
           #########   MOD SECTION FOR MANUAL SELECT INFLECTION POINT
           #infl.point <- 2000
           #plot.BCRanks[[nm]] <- qplot(bc_rank$total, bc_rank$rank, geom= "line") +
           #  geom_vline(xintercept = metadata(bc_rank)$knee, color="blue", linetype=2) +
           #  geom_vline(xintercept = infl.point, color = "green", linetype=2) +
           #  annotate("text", y=1000, x=1.5 *c(metadata(bc_rank)$knee, infl.point),
           #           label=c("knee", "inflection"), color = c("blue","green")) +
           #  scale_x_log10() + scale_y_log10() +
           #  labs(y="Barcode rank", x = "Total UMI count",title = paste0("Sample: ", nm))
           
           #res_mat[[nm]] <- res_mat.raw[[nm]][,totalCounts > infl.point]
           #######   MOD SECTION FOR MANUAL SELECT INFLECTION POINT
           
           res_mat[[nm]] <- res_mat.raw[[nm]][,totalCounts > bc_rank@metadata$inflection]
         }
         
         
  )
  ## read in transcriptome sparse matrix as Seurat object
  seu.list[[nm]] = CreateSeuratObject(counts = res_mat[[nm]], project = proj.name, min.cells = 3)
  #add sample name as metadata, later add other metadata values
  seu.list[[nm]]@meta.data$Sample <- nm
  
}  
#which sequencing technology was used
switch(seq_tech,
       "10x" = {
         pdf(paste(output.path.plots, "DropletUtils_isCells.pdf", sep = "/" ),height=5, width=5)
         for(nm in names(res_mat.raw) ){
           print(plot.BCRanks[[nm]])
         }
         dev.off()
       },
       
       "BD"= { 
         pdf(paste(output.path.plots, "BarcodeRanks_isCells.pdf", sep = "/" ),height=5, width=5)
         for(nm in names(res_mat.raw) ){
           print(plot.BCRanks[[nm]])
         }
         dev.off()
       },
       
       
       "Singleron"={
         pdf(paste(output.path.plots, "BarcodeRanks_isCells.pdf", sep = "/" ),height=5, width=5)
         for(nm in names(res_mat.raw) ){
           print(plot.BCRanks[[nm]])
         }
         dev.off()
       }
       
       
)
#rm(res_mat.raw,res_mat.raw)
###################   READ COUNT MATRIX   ######################
#Save raw seurat list
obj.name <- paste0("raw.",proj.name,".list.",today,".rds")
saveRDS(seu.list, file=paste(output.path,obj.name,sep="/") , 
        ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL)
#seu.list <- readRDS("/Volumes/MarkIV/178_sc_Gebaur/03_QC/11.04.2024/raw.scFFPE_Gebaur.list.11.04.2024.rds")
#plot the number of cells per sample
genes.boxplot <- data.frame() 
for(i in 1:length(seu.list)){
  aux <- seu.list[[i]][[c("Sample", "nFeature_RNA","nCount_RNA")]] %>% data.frame()
  genes.boxplot <- rbind(genes.boxplot,aux)
}
#Number of cells after kneeplot
df.cells.afterKnee <- genes.boxplot %>%
  group_by(Sample) %>%
  summarise(n = n())
###################       FIND DOUBLETS       ##################################
#Create doublet folder
output.path.scDblFinder <- paste(output.path,"scDblFinder",sep = "/")
dir.create(output.path.scDblFinder,recursive = TRUE)
seu.list.noDB <- seu.list
for(nm in names(seu.list.noDB)){
  #Transform Seurat object to SingleCellExperiment
  #sce <- as.SingleCellExperiment(seu.list.noDB[[nm]])
  #extract the counts as array of counts
  counts <- seu.list.noDB[[nm]]@assays$RNA@layers$counts
  colnames(counts) <- colnames(seu.list.noDB[[nm]])
  
  #Specifying doublet rate and standard deviation per capturing technology
  #which sequencing technology was used
  #switch(seq_tech,
  #       "10x" = {
  #         dbr=0.01
  #         dbr.sd= 0.015
  #         },
  #       "BD"= {#Need to create it
  #         dbr=NULL
  #         dbr.sd= NULL
  #         },
  #       "Singleron"={#Need to create it
  #         dbr=NULL
  #         dbr.sd= NULL
  #         }
  #)
  #Find doublets and tag them 
  #set.seed(1235)#For reproducible results
  #scDblFinder_output <- scDblFinder(counts,dbr=dbr, dbr.sd = dbr.sd, returnType= "table")
  scDblFinder_output <- scDblFinder(counts, returnType= "table")
  
  scDblFinder_output$BC <- rownames(scDblFinder_output)
  
  #Remove the artificial doublets by subseting only the ones tagged as "real"
  scDblFinder_output_real <- subset(scDblFinder_output, subset= type== "real")
  
  #Save results of Doublets to excel for future analysis
  writexl::write_xlsx(as.data.frame(scDblFinder_output), paste(output.path.scDblFinder, paste0("DoubletOuput_Sample.",nm,".xlsx") ,sep = "/"),col_names = TRUE)
  
  #Save possible doublet info to seurat
  seu.list.noDB[[nm]]$possDoublet <- scDblFinder_output_real$class
}
###################       FIND DOUBLETS       ##################################

#############       SUBSET SEURAT ONLY NO DOUBLET       ########################
for(nm in names(seu.list.noDB)){
  message(sprintf("Subsetting NO doublet sample= %s",nm))
  if ( "possDoublet" %in% names(seu.list.noDB[[nm]]@meta.data)) {
    #Filter only non doublet cells
    seu.list.noDB[[nm]] <- subset(seu.list.noDB[[nm]], subset = possDoublet == "singlet")
  }else{
    cat("Doublet finder did not run...")
    seu.list.noDB[[nm]]   <- seu.list[[nm]]
    next
  }
}
#############       SUBSET SEURAT ONLY NO DOUBLET       ########################
#plot the number of cells per sample
genes.boxplot <- data.frame() 
for(i in 1:length(seu.list.noDB)){
  aux <- seu.list.noDB[[i]][[c("Sample", "nFeature_RNA","nCount_RNA")]] %>% data.frame()
  genes.boxplot <- rbind(genes.boxplot,aux)
}
#Number of cells after kneeplot
df.cells.noDB <- genes.boxplot %>%
  group_by(Sample) %>%
  summarise(n = n())
#This seurat object is a raw seurat object after empty droplets and doublet removal. 
#Hence we must perfom the normal workflow for the QC
obj.name <- paste0("noDub.raw.",proj.name,".list.",today,".rds")
saveRDS(seu.list.noDB, file=paste(output.path,obj.name,sep="/") , 
        ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL)
#seu.list.noDB <- readRDS("/Volumes/MarkIV/178_sc_Gebaur/03_QC/11.04.2024/noDub.raw.scFFPE_Gebaur.list.11.04.2024.rds")


###################       NORMALIZING DATA       ###############################
seu.list.noDB.norm <- seu.list.noDB
for(i in 1:length(seu.list.noDB.norm)){
  if(specie=="hsa"){
    pat <- "^MT-"
  }else if(specie=="mmu"){
    pat <- "^mt-"
  }
  seu.list.noDB.norm[[i]][["percent.mt"]] <- PercentageFeatureSet(seu.list.noDB.norm[[i]],pattern = pat)
  
  #seu.list.norm[[i]] <- subset(seu.list.norm[[i]], subset = percent.mt < 25)
  seu.list.noDB.norm[[i]] <- NormalizeData(object = seu.list.noDB.norm[[i]],
                                           normalization.method = "LogNormalize",scale.factor = 10000)
  seu.list.noDB.norm[[i]] <- FindVariableFeatures(object = seu.list.noDB.norm[[i]],
                                                  selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  seu.list.noDB.norm[[i]] <- ScaleData(seu.list.noDB.norm[[i]], verbose=TRUE, rownames(seu.list.noDB.norm[[i]]))
  seu.list.noDB.norm[[i]] <- RunPCA(seu.list.noDB.norm[[i]], npcs=30, verbose=TRUE)
  if((seu.list.noDB.norm[[i]] %>% dim %>% '['(.,2)) > 150){
    seu.list.noDB.norm[[i]] <- RunTSNE(seu.list.noDB.norm[[i]], reduction = "pca", dims = 1:30, check_duplicates = FALSE)
    seu.list.noDB.norm[[i]] <- RunUMAP(seu.list.noDB.norm[[i]], reduction = "pca", dims = 1:30)
    seu.list.noDB.norm[[i]] <- FindNeighbors(seu.list.noDB.norm[[i]], reduction = "pca", dims = 1:30)
    seu.list.noDB.norm[[i]] <- FindClusters(seu.list.noDB.norm[[i]], resolution = 0.1, algorithm = 1, verbose = TRUE)
    
  } else {
    seu.list.noDB.norm[[i]] <- RunTSNE(seu.list.noDB.norm[[i]], reduction = "pca", dims = 1:30, check_duplicates = FALSE,perplexity =5)
    seu.list.noDB.norm[[i]] <- RunUMAP(seu.list.noDB.norm[[i]], reduction = "pca", dims = 1:30)
    seu.list.noDB.norm[[i]] <- FindNeighbors(seu.list.noDB.norm[[i]], reduction = "pca", dims = 1:30)
    seu.list.noDB.norm[[i]] <- FindClusters(seu.list.noDB.norm[[i]], resolution = 0.5, algorithm = 1, verbose = TRUE)
  }
}
###################       NORMALIZING DATA       ###############################
obj.name <- paste0("noDub.norm.",proj.name,".list.",today,".rds")
saveRDS(seu.list.noDB.norm, file=paste(output.path,obj.name,sep="/") , 
        ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL)
###################       VISUALIZATION   NO DB    ##################################
seu_plots(seu.obj.list=seu.list.noDB.norm, output.path =output.path.plots, processing.stage="02_NORM.NoDub")
###################       VISUALIZATION   NO DB    ##################################
#seu.list.noDB.norm <- readRDS("/Volumes/MarkIV/178_sc_Gebaur/03_QC/11.04.2024/noDub.norm.scFFPE_Gebaur.list.11.04.2024.rds")


######################      IDENTIFY LOW QUALITY CELLS  METHOD 2      ##################################
#THIS IS A MANUAL METHOD DONE BY A EXCEL FILE WHERE THE USER MUST SELECT THE PROPER THRESHOLDS FOR SUBSETING THE DATA
#Read QC metadata from file
sprintf("Reading QC filtering parameters from metadata file: %s", metadata_file)
QC_parameters <- readxl::read_excel(paste(orig,metadata_file,sep="/"),sheet="QC")
#Go through all of the samples
seu.list.noDB.norm.qc <- seu.list.noDB.norm
plots <- list()
for (nm in names(seu.list.noDB.norm.qc)){
  message(sprintf("Sample: %s", nm))
  #Extracting only the metadata from the specific sample
  Sample_QC_parameter <- QC_parameters %>% filter(Sample == nm)
  if(dim(Sample_QC_parameter)[1] != 0){ #check if metadata parameters were found
    #Create plots with the QC threshold plotted
    message(sprintf("Creating new plots with thresholds..."))
    plots[[nm]] <- seu_plots_withThresholdLines(seu.obj = seu.list.noDB.norm.qc[[nm]],qc_parameter_threshold = Sample_QC_parameter)
    
    #Subseting sample
    message(sprintf("Subseting sample with setted parameters..."))
    seu.list.noDB.norm.qc[[nm]] <- subset(seu.list.noDB.norm.qc[[nm]], subset= percent.mt < Sample_QC_parameter$MT_High )
    seu.list.noDB.norm.qc[[nm]] <- subset(seu.list.noDB.norm.qc[[nm]], subset= nCount_RNA > Sample_QC_parameter$nCount_RNA_Low & nCount_RNA < Sample_QC_parameter$nCount_RNA_High )
    seu.list.noDB.norm.qc[[nm]] <- subset(seu.list.noDB.norm.qc[[nm]], subset= nFeature_RNA > Sample_QC_parameter$nFeature_RNA_Low  & nFeature_RNA < Sample_QC_parameter$nFeature_RNA_High )
    #Saving QC parameters to misc
    seu.list.noDB.norm.qc[[nm]]@misc <- list("QC_Parameters"=list("MT"=paste0("X < ",Sample_QC_parameter$MT_High),
                                                                  "nCount_RNA"= paste(Sample_QC_parameter$nCount_RNA_Low,  Sample_QC_parameter$nCount_RNA_High, sep = " < X < "),
                                                                  "nFeature_RNA"=paste(Sample_QC_parameter$nFeature_RNA_Low, Sample_QC_parameter$nFeature_RNA_High, sep = " < X < ")
    )
    )
  }else {message("No QC metadata found")}
}
name <- paste0("03_NORM.NoDub_Pre-processed_Manual", "_ViolinPlot_Number_Cell.MT.Features_allSamples.pdf")
pdf(paste(output.path.plots, name , sep = "/" ))
print(plots)
dev.off()
######################      IDENTIFY LOW QUALITY CELLS  METHOD 2      ##################################
obj.name <- paste0("noDub.pre-processed.",proj.name,".list.",today,".rds")
saveRDS(seu.list.noDB.norm.qc, file=paste(output.path,obj.name,sep="/") , 
        ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL)
#plot the number of cells per sample
genes.boxplot <- data.frame() 
for(i in 1:length(seu.list.noDB.norm.qc)){
  aux <- seu.list.noDB.norm.qc[[i]][[c("Sample", "nFeature_RNA","nCount_RNA")]] %>% data.frame()
  genes.boxplot <- rbind(genes.boxplot,aux)
}
#Number of cells after kneeplot
df.cells.noDB.afterQC <- genes.boxplot %>%
  group_by(Sample) %>%
  summarise(n = n())


#Table number of cells
df.number.cells.list <- list(df.cells.afterKnee, df.cells.noDB, df.cells.noDB.afterQC)
df.number.cells <- df.number.cells.list %>% reduce(full_join, by='Sample')
colnames(df.number.cells) <- c("Sample","AfterKneePlot", "After Removing Multiplets", "AfterQC")
#pdf(paste(output.path.plots, paste0("Number_of_Cells.pdf"), sep = "/" ))
#grid.table(df.number.cells)
#dev.off()
#calculate total number of things
tot_samples <- paste0( length(df.number.cells$Sample) , " Samples")
tot_afterKnee <- sum(df.number.cells$AfterKneePlot) 
temp <- sum(df.number.cells$`After Removing Multiplets`)
tot_afterMultiplets <- paste0( temp, 
                               " ( ", 
                               round( (temp/tot_afterKnee)*100,digits = 2),
                               " % )")
temp2 <- sum(df.number.cells$AfterQC)
tot_afterQC <- paste0( temp2, 
                       " ( ", 
                       round( (temp2/tot_afterKnee)*100,digits = 2),
                       " % )")
tot <- c(tot_samples, tot_afterKnee,tot_afterMultiplets, tot_afterQC)
df.number.cells2 <- rbind(df.number.cells, tot)
#saving as pdf and xlsx
pdf(paste(output.path.plots, paste0("Number_of_Cells.pdf"), sep = "/" ))
grid.table(df.number.cells2)
dev.off()
#Saving number of cells in excel
#Filter_of_Viable_Cells <- c("Number of percent.mt < ", "Number of Features > ", "Number of Counts > ")
#value <- parameters_preprocessed
#df.parameters <- data.frame(Filter_of_Viable_Cells, value)
df.parameters <- data.frame()
for (nm in names(seu.list.noDB.norm.qc)){
  df.parameters <- rbind(df.parameters, cbind("Sample"= nm, as.data.frame(seu.list.noDB.norm.qc[[nm]]@misc) ))
}
wb <- createWorkbook()
## Add a worksheet
addWorksheet(wb, "Sheet 1", gridLines = FALSE)
writeData(wb, "Sheet 1", df.number.cells2, startCol = 1, startRow = 1, rowNames = FALSE, colNames = TRUE)
writeData(wb, "Sheet 1", df.parameters, startCol = 9, startRow = 1, rowNames = FALSE, colNames = TRUE)
saveWorkbook(wb,paste(output.path.plots, paste0("Number_of_Cells.xlsx"), sep = "/" ), overwrite = TRUE)

