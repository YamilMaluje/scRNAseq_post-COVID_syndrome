#This integration for SeuratV5
#This will read the processed Seurat list object from pre-processing_V2.X.R
#It will add the metadata to the cells by importing it from a xlsx file
#The user must set the metadata variable from which should be integrated
#It will integrate with different methods such as "cca","harmony", "rpca" for the user to compare which one you like best
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
options(future.globals.maxSize = 8000 * 1024^2)#o set it to 8GB
################# SETTING VARIABLES ##############################################
orig <- getwd()
srcRaw <- paste(orig,"03_QC","28.02.2025",sep = "/")
today <- Sys.Date() %>% format(., format="%d.%m.%Y")
#Create output folders to save all the files
output.path <- paste(orig,"04_integrated",today,sep = "/")
dir.create(output.path,recursive = TRUE)
proj.name <- "AFL_NEW_All"
specie<- "mmu" # homo sapiens "hsa" or mus musculus "mmu"
seu_processed_seuList <- "noDub.pre-processed.COVID-AFL_NEW_All.list.28.02.2025.rds"
metadata_file <- "metadata_covid_TNF_NEW_All.xlsx"
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
################     FUNCTIONS    ##############################################
integrate_different_methods <- function(obj,clu=0.5,split=NULL,methods="ALL",pca_thresh=100){
  #obj= is the merged object to be integrated
  #clu: resolution of cluster for FindClusters
  #split: Variable on the metadata used for spliting to fuurther remove the batch effect
  #methods: ALL if you want to run all of the supported integration methods, else write "cca","harmony", "rpca"
  library(future)
  method_names <- list("cca"="CCAIntegration",
                       "harmony"="HarmonyIntegration",
                       "rpca"="RPCAIntegration")
  
  if(methods=="ALL" || methods=="all" || methods=="A" || methods=="a"){
    cat("Running all integration methods...")
    for (nm in names(method_names)){
      print(nm)
      #Integrate with different methods
      obj <- IntegrateLayers(
        object = obj, 
        method =method_names[[nm]] ,
        orig.reduction = "pca", new.reduction = paste("integrated",nm,split,sep = "."),
        verbose = FALSE)
      obj <- FindNeighbors(obj, dims = 1:pca_thresh, reduction = paste("integrated",nm,split,sep = "."))
      obj <- FindClusters(obj, resolution = clu, algorithm = 4, verbose = TRUE, cluster.name = paste("integrated","clusters","res",clu,nm,split,sep = "."))
      obj <- RunUMAP(obj, dims = 1:pca_thresh, reduction = paste("integrated",nm,split,sep = "."), reduction.name = paste("umap","integrated",nm,split,sep = "."))
    }
  }else{
    
    cat("Running integration method: ", methods)
    obj <- IntegrateLayers(
      object = obj, 
      method =method_names[[methods]],
      orig.reduction = "pca", new.reduction = paste("integrated",methods,split,sep = "."),
      verbose = FALSE)
    
    obj <- FindNeighbors(obj, dims = 1:pca_thresh, reduction = paste("integrated",methods,split,sep = "."))
    obj <- FindClusters(obj, resolution = clu, algorithm = 4, verbose = TRUE, cluster.name = paste("integrated","clusters","res",clu,methods,split,sep = "."))
    obj <- RunUMAP(obj, dims = 1:pca_thresh, reduction = paste("integrated",methods,split,sep = "."), reduction.name = paste("umap","integrated",methods,split,sep = "."))
  }
  
  return(obj)
}
update_seurat_object_list <- function(obj.list.previous=NULL){
  #UPDATING FROM SEURAT V4 TO V5
  obj.list.previous.updated <- list()
  for (nm in names(obj.list.previous)){
    #obj.list.previous.updated[[nm]] <- UpdateSlots(object = obj.list.previous[[nm]])
    #obj.list.previous.updated[[nm]] <- UpdateSeuratObject(obj.list.previous.updated[[nm]])
    obj.list.previous.updated[[nm]] <- obj.list.previous[[nm]]
    # convert a v3 assay to a v5 assay
    obj.list.previous.updated[[nm]][["RNA"]] <- as(object = obj.list.previous[[nm]][["RNA"]], Class = "Assay5")
  }
  return(obj.list.previous.updated)
}
seurat_integration_by_variable <- function(merged_obj, 
                                           integration_variable="Sample",
                                           output_path=".",
                                           cluster_resolution=0.5, 
                                           methods="ALL"){
  #Create output folder for the plots and seurat object
  output.path_results <- paste(output_path,paste0("IntegratedBy_",integration_variable),sep = "/")
  dir.create(output.path_results,recursive = TRUE)
  message(sprintf("Output Folder created: %s", output.path_results))
  
  #Split the object by Sample
  merged_obj_split <- merged_obj
  DefaultAssay(merged_obj_split) <- "RNA"
  merged_obj_split@assays$SCT <- NULL
  aux <- unlist(merged_obj_split[[integration_variable]], use.names=TRUE) 
  names(aux) <- rownames(merged_obj_split[[integration_variable]])
  
  merged_obj_split[["RNA"]] <- split(merged_obj_split[["RNA"]], f = aux)
  #Layers(merged_obj_split[["RNA"]])
  
  #Processing the object for integration
  merged_obj_split <- NormalizeData(merged_obj_split)
  merged_obj_split <- FindVariableFeatures(merged_obj_split,nfeatures = 3000)
  merged_obj_split <- ScaleData(merged_obj_split)
  merged_obj_split <- RunPCA(merged_obj_split, npcs = 100, verbose = TRUE)
  merged_obj_split <- FindNeighbors(merged_obj_split, dims = 1:100, reduction = "pca")
  merged_obj_split <- FindClusters(merged_obj_split, resolution = cluster_resolution, algorithm = 4, verbose = TRUE, cluster.name = "unintegrated_clusters")
  merged_obj_split <- RunUMAP(merged_obj_split, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
  
  # visualize by Sample and cell type annotation
  p <- DimPlot(merged_obj_split, reduction = "umap.unintegrated", group.by = c( integration_variable ,"unintegrated_clusters"),label=T, repel = T) +
    labs(caption= paste0("Date of render: ",today))
  save_plot(path =output.path_results,
            plot_name = "00_UMAP_Unintegrated_BatchEffect",plot = p)
  
  
  # visualize by Sample and cell type annotation
  p <- DimPlot(merged_obj_split, reduction = "umap.unintegrated", split.by  = c(integration_variable),label=T, repel = T,raster = FALSE) +
    labs(caption= paste0("Date of render: ",today))
  save_plot(path =output.path_results,
            plot_name = "00_UMAP_Unintegrated_BatchEffect_split",plot = p)
  
  
  
  #merged_obj_split_original <- merged_obj_split
  #Integrate with different methods
  message(sprintf("Integrating seurat object by %s variable", integration_variable))
  integrated_obj_norm <- integrate_different_methods(obj=merged_obj_split,
                                                     clu=cluster_resolution,
                                                     split=integration_variable, 
                                                     methods = methods,
                                                     pca_thresh=30)
  
  
  
  # visualize by variable and check which integration method worked best
  p1 <- DimPlot(integrated_obj_norm, reduction = "umap.unintegrated", group.by= c(integration_variable),label=T, repel = T, raster = FALSE) +
    labs(caption= paste0("Date of render: ",today))
  
  temp <- names(integrated_obj_norm@reductions)
  # Use grep to find indices of matching elements
  matching_indices <- grep(paste("umap","integrated","*",integration_variable,sep = "."), temp)
  # Use the indices to subset the vector
  redu_names <- temp[matching_indices]
  p <- list()
  for (redu in redu_names){
    p[[redu]] <- DimPlot(integrated_obj_norm, reduction = redu, group.by= integration_variable, label=T, repel = T,raster = FALSE) +
      labs(caption= paste0("Date of render: ",today))
  }
  plot_name <- paste(output.path_results,"01_UMAP_Integrated_AllMethods.pdf", sep = "/")
  plots <- CombinePlots(
    plots = c(list(p1),p),
    ncol = 2)
  pdf(file= plot_name,onefile = TRUE, width = 12, height = 8)
  print(plots)
  dev.off()
  
  # re-join layers after integration
  integrated_obj_norm[["RNA"]] <- JoinLayers(integrated_obj_norm[["RNA"]])
  
  #Save the integration object
  integrated_file_path <- paste(output.path_results,paste0(proj.name,".integratedby",integration_variable,"_",today,".rds"),sep="/")
  sprintf("Saving integrated object in: %s",integrated_file_path)
  saveRDS(integrated_obj_norm, file=integrated_file_path)  
  
}
save_plot <- function(path=output.path, plot=NULL, plot_name=NULL){
  plot_name <- paste0(paste(plot_name,today,sep = "_"), ".pdf")
  ggsave(filename= plot_name , plot=plot, path=path, width = 14, height = 8, dpi = 300) 
}

################     FUNCTIONS    ##############################################
#################    READ SEURAT LIST     ######################################
#Read seurat object list all samples
obj.list <- readRDS(paste(srcRaw,seu_processed_seuList,sep = "/"))
#UNCOMMENT IF YOU WANT TO ADD OBJECT LIST OF PREVIOUS SEURAT VERSIONS
#Read Previous object list as well to be added to the new samples
#obj.list.previous <- readRDS("/Users/yamilmac/Documents/Projects/UKSH scRNAseq/Server Files/scRNAseq Analysis/Köhl/sc_Köhl_mouse/Samyr_normlist_noDub_09.06.2023.rds")
format(object.size(obj.list), units = "Mb")
#################    READ SEURAT LIST     ######################################
#UNCOMMENT IF YOU WANT TO ADD OBJECT LIST OF PREVIOUS SEURAT VERSIONS
#################    APPENDING SUERAT OBJECT     ###############################
#obj.list.previous.updated <- update_seurat_object_list(obj.list.previous = obj.list.previous)
#new_obj_list <- append(obj.list, obj.list.previous.updated)
#################    APPENDING SUERAT OBJECT     ###############################
################# ADDING METADATA TO SEURAT LIST  ##############################
#READ metadata from file
sprintf("Reading QC filtering parameters from metadata file: %s", metadata_file)
new.meta <- readxl::read_excel(paste(orig,metadata_file,sep="/"),sheet="Meta")
for (nm in names(obj.list)){
  message("Adding metadata to sample ", nm)
  #Extracting only the metadata from the specific sample
  sample.meta <- new.meta %>% filter(Sample == nm)
  if(dim(sample.meta)[1] != 0){
    #Adding the metadata to the seurat list with the same name as the column name
    for (col.name in colnames(sample.meta)) {
      #message("Adding ", col.name, "= ",sample.meta[[col.name]] , " to sample ", nm)
      obj.list[[nm]][[col.name]] <- sample.meta[[col.name]]
    }
  }else {message("No metadata found")}
}
################# ADDING METADATA TO SEURAT LIST  ##############################
################# MERGING OBJECTS WITHOUT INTEGRATION  #########################
#Merging the Seurat list samples
aux <- obj.list[-1] #will remove the first seurat object 
cell_id <- names(obj.list) #NAmes of the samples for future index id
merged_obj <- merge(obj.list[[1]], y = aux , add.cell.ids = cell_id, project = proj.name,  merge.data = TRUE)
#If not proceeding with integration, rejoin the layers after merging.
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])
Layers(merged_obj[["RNA"]])
rm(aux)
#Save the merged object
merged_file_path <- paste(output.path,paste0(proj.name,".merged_od_norm.",today,".rds"),sep="/")
sprintf("Saving integrated object in: %s",merged_file_path)
saveRDS(merged_obj, file=merged_file_path)
################# MERGING OBJECTS WITHOUT INTEGRATION  #########################


################################################################################
################# INTEGRATING BY DIFFERENT VARIABLES  ##########################
################################################################################
#Subsample for testing 
# Downsample the number of cells per identity class
#Idents(merged_obj) <- "Sample"
# Downsample the number of cells per identity class
#merged_obj_sub <- subset(x = merged_obj, downsample = 100)
#Sketch dataset to test a representative sample out of the full dataset
#https://satijalab.org/seurat/articles/seurat5_sketch_analysis
#Integrate by Batch
seurat_integration_by_variable(merged_obj, 
                               integration_variable="Batch",
                               output_path=output.path, 
                               cluster_resolution=0.1, 
                               methods="harmony")
################################################################################
################# INTEGRATING BY DIFFERENT VARIABLES  ##########################
########################################################################### 
