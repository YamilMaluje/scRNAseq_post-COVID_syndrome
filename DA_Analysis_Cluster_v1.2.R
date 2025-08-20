################################################################################
################################################################################
#############           DASEQ subpopulation                #####################
################################################################################
################################################################################
#Paper: https://www.pnas.org/doi/10.1073/pnas.2100293118
#Website: https://github.com/KlugerLab/DAseq
#Tutorial: https://klugerlab.github.io/DAseq/articles/tutorial.html
#Install library
#devtools::install_github("KlugerLab/DAseq")
#Import Library
#library(DAseq)
library(dplyr)
library(ggplot2)  
library(stringr)
library(Seurat)
#library(EnhancedVolcano)
#Cluster
library(reticulate)
use_condaenv("DAseq")
python2use <- "/home/maluje/.local/share/r-miniconda/envs/DAseq/bin/python"
GPU <- 1
version <- 1.1
description <- "Analysis of Differentially Abundance in a Seurat object..."
sprintf("Version= %.1f  %s ",version, description)
library(DAseq)
################################################################################
################################################################################
##############################  INPUT ARGUMENTS ################################
args = commandArgs(TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args) >= 1) {
  
  #PAth to the list of pre-processed seurat object
  #Current working directory
  seu.obj.path = args[1]
  #Output directory to save all the analysis
  output.dir = args[2]
  #Variabble of comparison
  #i.e: Cancer:0:1
  #This example will subset the seurat obbject by "Cancer" and will assign the cancer cells that has "0" to the Blue 
  # and the cells that has "1" to the Red color
  variable.comparison = args[3]
  #reduction to use for plots (i.e. UMAP, TSNE), could be a string or a string vector
  reductions =args[4]
  
}
sprintf("INPUT ARGUMENTS:")
sprintf("seurat object path: %s",seu.obj.path)
sprintf("Absolute output path where the plots will be saved : %s",output.dir)
sprintf("Variable of Comparisonn : %s",variable.comparison)
sprintf("Dimensional reduction plots : %s",reductions)
##############################  INPUT ARGUMENTS ################################
################################################################################


################################################################################
##################           FUNCTIONS                     #####################
################################################################################
getDAregion<-function (X, da.cells, cell.labels, labels.1, labels.2, prune.SNN = 1/15, 
    resolution = 0.05, group.singletons = F, min.cell = NULL, 
    do.plot = T, plot.embedding = NULL, size = 0.5, do.label = F, 
    ...) 
{
    if (!inherits(x = X, what = "matrix")) {
        cat("Turning X to a matrix.\n")
        X <- as.matrix(X)
    }
    if (!inherits(cell.labels, "character") | !inherits(labels.1, 
        "character") | !inherits(labels.2, "character")) {
        stop("Input parameters cell.labels, labels.1 and labels.2 must be character")
    }
    if (length(setdiff(cell.labels, c(labels.1, labels.2))) > 
        0) {
        warning("Input parameter cell.labels contain labels not from labels.1 or labels.2, subsetting...")
        cell.idx <- which(cell.labels %in% c(labels.1, labels.2))
        X <- X[cell.idx, ]
        cell.labels <- cell.labels[cell.idx]
    }
    else {
        cell.idx <- seq_len(length(cell.labels))
    }
    if (is.null(min.cell)) {
        min.cell <- as.integer(colnames(da.cells$da.ratio)[1])
        cat("Using min.cell = ", min.cell, "\n", sep = "")
    }
    n.cells <- nrow(X)
    n.dims <- ncol(X)
    if (is.null(rownames(X))) {
        rownames(X) <- paste("C", c(1:n.cells), sep = "")
    }
    seurat.version <- substr(packageVersion("Seurat"), 1, 1)
    if (seurat.version == "3" | seurat.version == "4" | seurat.version == "5") { #UPDATED
        X.S <- CreateSeuratObject(counts = t(X))
        #X.S@reductions$pca <- new("DimReduc", cell.embeddings = X, 
        #    assay.used = DefaultAssay(X.S), key = "PC_")
        X.S@reductions$pca <- CreateDimReducObject(embeddings = X, #UPDATED
            assay = DefaultAssay(X.S), key = "PC_")
        X.S <- FindNeighbors(X.S, reduction = "pca", dims = 1:n.dims, 
            prune.SNN = prune.SNN, verbose = F)
        if (length(da.cells$da.up) > 1) {
            up.S <- CreateSeuratObject(counts = t(X[da.cells$da.up, 
                ]))
            #up.S@reductions$pca <- new("DimReduc", cell.embeddings = X[da.cells$da.up, 
            #    ], assay.used = DefaultAssay(up.S), key = "PC_")
            up.S@reductions$pca <- CreateDimReducObject(embeddings = X[da.cells$da.up, #UPDATED
                ], assay = DefaultAssay(up.S), key = "PC_")
            up.S <- FindNeighbors(up.S, reduction = "pca", dims = 1:n.dims, 
                verbose = F)
            up.snn <- X.S@graphs$RNA_snn[da.cells$da.up, da.cells$da.up]
            up.S@graphs$RNA_snn@i <- up.snn@i
            up.S@graphs$RNA_snn@p <- up.snn@p
            up.S@graphs$RNA_snn@x <- up.snn@x
            #up.S <- FindClusters(up.S, resolution = resolution, 
            #    group.singletons = group.singletons, verbose = F, ...)
            up.S <- FindClusters(up.S, resolution = resolution, #UPDATED
                verbose = F)
            up.clusters <- as.numeric(up.S@active.ident)
            up.clusters[up.S@active.ident == "singleton"] <- 0
        }
        else {
            up.clusters <- NULL
        }
        n.up.clusters <- length(unique(up.clusters)) - as.numeric(0 %in% 
            up.clusters)
        if (length(da.cells$da.down) > 1) {
            down.S <- CreateSeuratObject(counts = t(X[da.cells$da.down, 
                ]))
            #down.S@reductions$pca <- new("DimReduc", cell.embeddings = X[da.cells$da.down, 
            #    ], assay.used = DefaultAssay(down.S), key = "PC_")
            down.S@reductions$pca <- CreateDimReducObject(embeddings = X[da.cells$da.down, #UPDATE
                ], assay = DefaultAssay(down.S), key = "PC_")
            down.S <- FindNeighbors(down.S, reduction = "pca", 
                dims = 1:n.dims, verbose = F)
            down.snn <- X.S@graphs$RNA_snn[da.cells$da.down, 
                da.cells$da.down]
            down.S@graphs$RNA_snn@i <- down.snn@i
            down.S@graphs$RNA_snn@p <- down.snn@p
            down.S@graphs$RNA_snn@x <- down.snn@x
            #down.S <- FindClusters(down.S, resolution = resolution, 
            #    group.singletons = group.singletons, verbose = F, 
            #    ...)
            down.S <- FindClusters(down.S, resolution = resolution, #UPDATED
                verbose = F)
            down.clusters <- as.numeric(down.S@active.ident) + 
                n.up.clusters
            down.clusters[down.S@active.ident == "singleton"] <- 0
        }
        else {
            down.clusters <- NULL
        }
    }
    da.region.label <- rep(0, n.cells)
    da.region.label[da.cells$da.up] <- up.clusters
    da.region.label[da.cells$da.down] <- down.clusters
    da.region.label.tab <- table(da.region.label)
    if (min(da.region.label.tab) < min.cell) {
        da.region.to.remove <- as.numeric(names(da.region.label.tab)[which(da.region.label.tab < 
            min.cell)])
        cat("Removing ", length(da.region.to.remove), " DA regions with cells < ", 
            min.cell, ".\n", sep = "")
        da.region.label.old <- da.region.label
        for (ii in da.region.to.remove) {
            da.region.label[da.region.label.old == ii] <- 0
            da.region.label[da.region.label.old > ii] <- da.region.label[da.region.label.old > 
                ii] - 1
        }
    }
    X.n.da <- length(unique(da.region.label)) - 1
    X.da.stat <- matrix(0, nrow = X.n.da, ncol = 3)
    colnames(X.da.stat) <- c("DA.score", "pval.wilcoxon", "pval.ttest")
    if (X.n.da > 0) {
        for (ii in 1:X.n.da) {
            X.da.stat[ii, ] <- getDAscore(cell.labels = cell.labels, 
                cell.idx = which(da.region.label == ii), labels.1 = labels.1, 
                labels.2 = labels.2)
        }
    }
    else {
        warning("No DA regions found.")
    }
    if (do.plot & is.null(plot.embedding)) {
        warning("plot.embedding must be provided by user if do.plot = T")
        X.region.plot <- NULL
    }
    else if (do.plot & !is.null(plot.embedding)) {
        plot.embedding <- plot.embedding[cell.idx, ]
        X.da.label <- da.region.label
        X.da.order <- order(X.da.label, decreasing = F)
        X.region.plot <- plotCellLabel(X = plot.embedding[X.da.order, 
            ], label = as.factor(X.da.label[X.da.order]), size = size, 
            do.label = do.label, label.plot = as.character(c(1:X.n.da))) + 
            scale_color_manual(values = c("gray", scales::hue_pal()(X.n.da)), #UPDATED
                breaks = c(1:X.n.da))
    }
    else {
        X.region.plot <- NULL
    }
    return(list(cell.idx = cell.idx, da.region.label = da.region.label, 
        DA.stat = X.da.stat, da.region.plot = X.region.plot))
}
getDAscore <- function(cell.labels, cell.idx, labels.1, labels.2){
  labels.1 <- labels.1[labels.1 %in% cell.labels]
  labels.2 <- labels.2[labels.2 %in% cell.labels]
  idx.label <- cell.labels[cell.idx]
  ratio.1 <- sum(idx.label %in% labels.1) / sum(cell.labels %in% labels.1)
  ratio.2 <- sum(idx.label %in% labels.2) / sum(cell.labels %in% labels.2)
  ratio.diff <- (ratio.2 - ratio.1) / (ratio.2 + ratio.1)
  cell.label.name <- sort(unique(cell.labels))
  cell.label.tab <- table(factor(cell.labels, levels = cell.label.name))
  idx.label.ratio <- table(factor(idx.label, levels = cell.label.name)) / cell.label.tab
  # print(idx.label.ratio)
  # score <- (mean(idx.label.ratio[labels.2]) - mean(idx.label.ratio[labels.1]))
  # score.n <- (mean(idx.label.ratio[labels.2]) - mean(idx.label.ratio[labels.1])) / sum(idx.label.ratio)
  if(length(labels.1) > 1 & length(labels.2) > 1){
    pval.wilcox <- wilcox.test(x = idx.label.ratio[labels.2], idx.label.ratio[labels.1])$p.value
    pval.ttest <- t.test(x = idx.label.ratio[labels.2], idx.label.ratio[labels.1])$p.value
  } else {
    pval.wilcox <- NA
    pval.ttest <- NA
  }
  return(c("DA.score" = ratio.diff, "pval.wilcoxon" = pval.wilcox, "pval.ttest" = pval.ttest))
}
################################################################################
##################           FUNCTIONS                     #####################
################################################################################


################################################################################
##################           MAIN                     ##########################
################################################################################
variable.comparison.vector <- str_split_fixed(variable.comparison, ":", 3) %>% as.vector(.)
sprintf("Comparison of %s : %s Blue and %s Red",variable.comparison.vector[1],variable.comparison.vector[2],variable.comparison.vector[3] )
sprintf("Load seurat object...%s",seu.obj.path)
obj <- readRDS(seu.obj.path)
sprintf("Subset only the cells that are in the %s",variable.comparison.vector[1])
#Subset only the cells that are in the variable.comparison.vector[1]
expr <- FetchData(obj, vars = variable.comparison.vector[1])
obj.sub <- obj[, which(x=   expr == variable.comparison.vector[2] |
                         expr == variable.comparison.vector[3] ) ]
sprintf("Retrieve sample names from the subset dataset")
#sample label names that correspond to one biological condition 
label_no <- obj[, which(x=   expr == variable.comparison.vector[2] )] %>% .$Sample %>%  unique(.)
label_yes <-obj[, which(x=   expr == variable.comparison.vector[3] )] %>% .$Sample %>%  unique(.)
#All sample labels of all cells
labels.all <- obj.sub$Sample %>% as.vector(.)

for (reduction in reductions){
  message(sprintf("Running the DA with %s reduction",reduction))
  
  #Creatinng the output directory
  output.dir.reduction <- paste(output.dir,variable.comparison.vector[1],reduction,sep = "/")
  message(sprintf("Creatinng the output directory %s",output.dir.reduction))
  dir.create(output.dir.reduction,recursive = TRUE)
  
  message(sprintf("COMPUTE DIFFERENTIALLY ABUNDANCE (DA) ON ALL CELLS"))
  # 01  ######  COMPUTE DIFFERENTIALLY ABUNDANCE (DA) ON ALL CELLS #######
  da_cells <- getDAcells(
    X = obj.sub@reductions[["pca"]]@cell.embeddings[,1:10],
    cell.labels = labels.all,
    labels.1 = label_no, #Blue | da.down
    labels.2 = label_yes, #Red | da.up
    k.vector = seq(50, 500, 50),
    plot.embedding = obj.sub@reductions[[reduction]]@cell.embeddings
  )
  
  #The prediction values are overlayed on the 2D embedding in the pred.plot slot of the output.
  plot <- da_cells$pred.plot + 
    ggtitle(paste0(reduction, " (Prediction Values)"))
  ggsave(filename= paste0(reduction,"_prediction_value.pdf") , plot=plot, path=output.dir.reduction) 
  
  #Random permutation on the labels to generate the null distribution for DA measure. 
  #We can check the permutation results in the rand.plot slot.
  plot <- da_cells$rand.plot +
    ggtitle("Permutation Results", subtitle = "Random permutation on the labels to generate the null distribution for DA measure")
  ggsave(filename= "Permutation_results.pdf" , plot=plot, path=output.dir.reduction) 
  
  #Selected DA cells are highlighted in the 2D embedding tsne
  #BLUE : these cells are more abundant in variable.comparison.vector[1] == variable.comparison.vector[2]
  #RED :  these cells are more abundant in variable.comparison.vector[1] == variable.comparison.vector[3]
  plot <- da_cells$da.cells.plot +
    ggtitle(paste0("Selected differentially abundant cells [", variable.comparison.vector[1] ,"]"), 
            subtitle = paste0("Red: ",variable.comparison.vector[3]," | Blue: ",variable.comparison.vector[2])
    )
  ggsave(filename= "DA_cells.pdf" , plot=plot, path=output.dir.reduction,width = 12, height = 10, dpi = 200 ) 
  
  #Save RDS and load it
  saveRDS(da_cells, paste(output.dir.reduction,paste0("da_cells_",reduction,".rds"),sep = "/"))
  
  
  message(sprintf("COMPUTE DIFFERENTIALLY ABUNDANCE REGIONS"))
  # 02 ######  COMPUTE DIFFERENTIALLY ABUNDANCE REGIONS #######
  #Selected DA cells will be clustered into several coherent regions, which represent the DA cell subpopulations.
  da_regions <- getDAregion(
    X = obj.sub@reductions[["pca"]]@cell.embeddings[,1:10],
    da.cells = da_cells,
    cell.labels = labels.all,
    labels.1 = label_no, #Blue | da.down
    labels.2 = label_yes, #Red | da.up
    resolution = 0.05,
    plot.embedding = obj.sub@reductions[[reduction]]@cell.embeddings,
  )
  
  #Clustering result is shown in the plot: da.region.plot slot
  plot<- da_regions$da.region.plot +
    ggtitle("Differentially abundant clusters", subtitle = "Represent the DA cell subpopulations")
  ggsave(filename= paste0("DA_",reduction,"_cluster.pdf") , plot=plot, path=output.dir.reduction, width = 12, height = 10, dpi = 200 ) 
  
  #Save RDS and load it
  saveRDS(da_regions, paste(output.dir.reduction,paste0("da_regions_",reduction,".rds"),sep = "/"))
  
  message(sprintf("Get markers for each DA subpopulation with STG"))
  # 03 ####    Get markers for each DA subpopulation with STG #####
  #characterize each DA subpopulation by detecting genes that seprate the DA subpopulation from the rest of the cells through STG (stochastic gates).
  
  #Stochastic gates (STG) to identify markers for each DA subpopulation.
  # select genes that separate each DA region from the rest of the cells
  STG_markers <- STGmarkerFinder(
    X = obj.sub[["RNA"]]$data,
    da.regions = da_regions,
    lambda = 1.5, n.runs = 5, return.model = T,
    python.use = python2use, GPU = GPU
  )
  #Save RDS and load it
  saveRDS(STG_markers, paste(output.dir.reduction,paste0("STG_markers_",reduction,".rds"),sep = "/"))
  
  #plot the prediction value from STG for DA subpopulation 
  plot.list <- list()
  for (i in 1:length(STG_markers$da.markers)){
    
    plot.list[[i]]<- plotCellScore(
      X = obj.sub@reductions[[reduction]]@cell.embeddings,
      score = STG_markers$model[[i]]$pred) +
      ggtitle("Markers score for each DA subpopulation", subtitle = paste0("Subpopulations: ",i))
    
    name <- paste0("Marker_DA_subpopulation_",i,".pdf")
    ggsave(filename= name , plot=plot.list[[i]], path=output.dir.reduction) 
  }
  
  
}
################################################################################
##################           MAIN                     ##########################
################################################################################
