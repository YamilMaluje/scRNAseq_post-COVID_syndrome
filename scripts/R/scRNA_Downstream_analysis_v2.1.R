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
################# SETTING VARIABLES ##############################################
orig <- getwd()
srcRaw <- paste(orig,"04_integrated","09.04.2025","IntegratedBy_Batch",sep = "/")
today <- Sys.Date() %>% format(., format="%d.%m.%Y")
output.path <- paste(orig,"05_DownstreamAnalysis",today,sep = "/")
output.path.plots <- paste(output.path, "Plots",sep = "/")
dir.create(output.path.plots,recursive = TRUE)
clu <- c(0.1)
proj.name <- "AFL_NEW_All"
specie<- "hsa" # homo sapiens "hsa" or mus musculus "mmu"
integrated_by <- "Batch"
integrated_file_name <- "AFL_NEW_All.integratedbyBatch_09.04.2025.rds" 
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
#######################    FUNCTIONS      ######################################

#################    READ SEURAT LIST     ######################################
#Read seurat object list all samples
integrated <- readRDS(paste(srcRaw, integrated_file_name,sep="/"))
#################    READ SEURAT LIST     ######################################



###################            VIZUALIZE DATA           ########################
output.path.plots.general <- paste(output.path.plots, "01_General",sep = "/")
dir.create(output.path.plots.general,recursive = TRUE)
p <- DimPlot(integrated,
             reduction = "umap.integrated.harmony.Batch",
             group.by = "Batch",
             cols=custom_colors,raster = FALSE) + 
  ggtitle(label = "UMAP group by Batch") 
save_plot(path = output.path.plots.general, plot_name = "UmapClusteringBatch",plot = p)
p <- DimPlot(integrated,
             reduction = "umap.integrated.harmony.Batch",
             split.by = "Batch",
             cols=custom_colors,raster = FALSE) + 
  ggtitle(label = "UMAP split by Batch") 
save_plot(path = output.path.plots.general, plot_name = "UmapClusteringBatch_split",plot = p)

p <- DimPlot(integrated,
             reduction = "umap.integrated.harmony.Batch",
             split.by = "Sample",
             cols=custom_colors,raster = FALSE) + 
  ggtitle(label = "UMAP split by Sample") 
save_plot(path = output.path.plots.general, plot_name = "UmapClusteringSample_split",plot = p)

p = DimPlot(integrated,
            reduction = "umap.integrated.harmony.Batch",
            split.by = "Condition",
            label=TRUE,
            cols=custom_colors,
            repel = TRUE,
            raster = FALSE) + 
  ggtitle(label = "UMAP split by Condition") 
save_plot(path = output.path.plots.general, plot_name = "UmapClusteringConditionSplit",plot = p)


#Plot number of cells per cluster
name.ident <- paste("integrated.clusters.res",clu,"harmony", integrated_by,sep=".")
dat <- aggregate(list("Frequency" = integrated@meta.data[[name.ident]]), 
                 integrated@meta.data[c(name.ident)],
                 FUN = length)
colnames(dat) <- c("X","Frequency")
rownames(dat) <- dat$X
dat <- dat[order(as.numeric(rownames(dat))), , drop = FALSE]



g<- ggplot(dat, aes( x=reorder(X,-Frequency), y=Frequency))
g  <- g + geom_col()+
  geom_text(aes(label = Frequency, vjust= -0.5),size=3)+
  labs(title="Distribution of cells per cell cluster", 
       subtitle=paste0("Resolution: ",clu),
       x="Cluster",
       y="Cells",
       caption= paste0("Date of render: ",today))
save_plot(path = output.path.plots.general, plot_name = paste0("CellDistributionPerCluster_res",clu), plot = g)



#Plot distribution of sample within each cell cluster
p2.data <- list()
library(stringr)
for (res in clu){ 
  dat <- aggregate(list("Frequency" = integrated@meta.data[["Sample"]]), integrated@meta.data[c("Sample",name.ident)], FUN = length)
  
  #split the sample name to further order it
  dat$Sample <- paste("Sample",dat$Sample,sep = "_")
  colnames(dat) <- c("Sample","X","Frequency")
  p2.data[[paste0("res_",res)]] <- ggplot(dat,aes(fill=Sample, y=Frequency, x=X)) +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual("legend", values = custom_colors) +
    labs(title="Distribution of sample within each cell cluster", 
         subtitle=paste0("Resolution: ",res),
         x="Cluster",
         y="Sample",
         caption= paste0("Date of render: ",today)
    )
  save_plot(path = output.path.plots.general, plot_name = paste0("SampleDistributionPerCluster_res",res), plot = p2.data[[paste0("res_",res)]])
  
}

p2.data <- list()
library(stringr)
for (res in clu){ 
  dat <- aggregate(list("Frequency" = integrated@meta.data[["Condition"]]), integrated@meta.data[c("Condition",name.ident)], FUN = length)
  
  #split the sample name to further order it
  dat$Condition <- paste("Condition",dat$Condition,sep = "_")
  colnames(dat) <- c("Condition","X","Frequency")
  p2.data[[paste0("res_",res)]] <- ggplot(dat,aes(fill=Condition, y=Frequency, x=X)) +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual("legend", values = custom_colors) +
    labs(title="Distribution of sample within each cell cluster", 
         subtitle=paste0("Resolution: ",res),
         x="Cluster",
         y="Condition",
         caption= paste0("Date of render: ",today)
    )
  save_plot(path = output.path.plots.general, plot_name = paste0("ConditionDistributionPerCluster_res",res), plot = p2.data[[paste0("res_",res)]])
  
}


p2.data <- list()
library(stringr)
for (res in clu){ 
  #name.ident <- "integrated.clusters.harmony.Batch"
  dat <- aggregate(list("Frequency" = integrated@meta.data[["Batch"]]), integrated@meta.data[c("Batch",name.ident)], FUN = length)
  
  #split the sample name to further order it
  dat$Batch <- paste("Batch",dat$Batch,sep = "_")
  colnames(dat) <- c("Batch","X","Frequency")
  p2.data[[paste0("res_",res)]] <- ggplot(dat,aes(fill=Batch, y=Frequency, x=X)) +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual("legend", values = custom_colors) +
    labs(title="Distribution of sample within each cell cluster", 
         subtitle=paste0("Resolution: ",res),
         x="Cluster",
         y="Batch",
         caption= paste0("Date of render: ",today)
    )
  save_plot(path = output.path.plots.general, plot_name = paste0("BatchDistributionPerCluster_res",res), plot = p2.data[[paste0("res_",res)]])
  
}
###################            VIZUALIZE DATA           ########################


########################      FIND MARKERS        ##############################
#Finding all marker for all cluster resolution using MAST and negbinomial Test
#clu <- c(0.5)
for (res in clu){
  message("Setting default idents res ",res, " and RNA assay...")
  Idents(integrated)<- name.ident
  DefaultAssay(integrated) <- "RNA"
  
  #Save the find markers in the misc slot of each object list
  #use MAST as the test 
  message("Finding all markers of res ", res, " using MAST test...")
  FindAllMarkersMAST <- FindAllMarkers(integrated,
                                       assay= "RNA",
                                       test.use="MAST",
                                       only.pos = FALSE, 
                                       min.pct = 0.25, 
                                       logfc.threshold = 0.25)
  slot.name <- paste0("FindAllMarkersMAST.res",res)
  message("Saving results in ",slot.name)
  integrated@misc[[slot.name]] <- FindAllMarkersMAST
  
  writexl::write_xlsx(FindAllMarkersMAST,path = paste(output.path.plots.general,
                                                      paste0(slot.name,".xlsx"),
                                                      sep = "/"))
  
  #message("Finding all markers of res ", res, " using negbinom test")
  #FindAllMarkersNegbinom <- FindAllMarkers(integrated,
  #                                         assay= "RNA",
  #                                         test.use="negbinom",
  #                                         only.pos = FALSE, 
  #                                         min.pct = 0.25, 
  #                                         logfc.threshold = 0.25)
  #slot.name <- paste0("FindAllMarkersNegbinom.res",res)
  #message("Saving results in ",slot.name, "\n")
  #integrated@misc[[slot.name]] <- FindAllMarkersNegbinom
  
  #writexl::write_xlsx(FindAllMarkersNegbinom,path = paste(output.path.plots.general,
  #                                                        paste0(slot.name,".xlsx"),
  #                                                        sep = "/"))
  
  
}
#add from information from which seurat processed object was integrated
integrated@misc["Findmarkers_from_"] <- integrated_file_name
#Save the integration object
integrated_file_path <- paste(output.path,paste0(proj.name,".integrated_od_markers_",today,".rds"),sep="/")
sprintf("Saving integrated object in: %s",integrated_file_path)
saveRDS(integrated, file=integrated_file_path)
########################      FIND MARKERS        ##############################



########################      FIND CONSERVED MARKERS        ####################
#Finding all marker for all cluster resolution using MAST and negbinomial Test
#clu <- c(0.5)
for (res in clu){
  message("Setting default idents res ",res, " and RNA assay...")
  Idents(integrated)<- name.ident
  DefaultAssay(integrated) <- "RNA"
  
  conserved_markers <- list()
  for(i in levels(integrated@meta.data[[name.ident]])  ){
    conserved_markers[[i]] <- FindConservedMarkers(integrated, ident.1 = i , grouping.var = "Condition", verbose = FALSE)
  }
  slot.name <- paste0("Conserved_markers.res",res)
  writexl::write_xlsx(conserved_markers,path = paste(output.path.plots.general,
                                                     paste0(slot.name,".xlsx"),
                                                     sep = "/"))
  
}
#Setting default idents res 0.1 and RNA assay...
#Warning: Identity: 2 not present in group TNFa_TGFb. Skipping TNFa_TGFb
#Warning: Identity: 3 not present in group TNFa_TGFb. Skipping TNFa_TGFb
#Warning: Identity: 3 not present in group TGFb. Skipping TGFb
#Warning: TNFa_TGFb has fewer than 3 cells in Identity: 5. Skipping TNFa_TGFb
#Warning: Identity: 6 not present in group PBS. Skipping PBS
########################      FIND CONSERVED MARKERS        ####################





########################     DGE ANALYSIS         ##############################
#integrated <- readRDS("")
output.path.plots.DGE <- paste(output.path, "Plots","02_DGE",sep = "/")
dir.create(output.path.plots.DGE,recursive = TRUE)
#For Differential gene expression (FindAllMarkers) "RNA" assay must be set to default
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)
integrated <- FindVariableFeatures(integrated)
integrated <- ScaleData(integrated)
#clu <- c("0.1","0.3","0.5", "0.8")
heatmap.list <- list()
topN <- 10
for (res in clu){
  message("Setting default idents res ",res, " and RNA assay...")
  Idents(integrated)<- name.ident
  Idents(integrated) <- factor(x = Idents(integrated), levels = sort(as.numeric(levels(integrated))))
  DefaultAssay(integrated) <- "RNA"
  
  # Select unique top10 genes resolution 0.1
  test1 <- paste0("FindAllMarkersMAST.res",res)
  #top <-  integrated@misc[[test1]]  %>% group_by(cluster) %>% top_n(topN, avg_log2FC)
  top <-  integrated@misc[[test1]]  %>% 
    group_by(cluster) %>% 
    arrange(.,  p_val_adj, desc(avg_log2FC))%>%
    slice_head(n = topN)
  
  #utop<-unique(top$gene)
  #Heatmap of the DGE
  heatmap.list[[test1]] <- DoHeatmap(integrated, features = top$gene) +    
    ggtitle(label = "DGE [MAST]",subtitle = paste0("Resolution: ",res)) +
    labs(caption= paste0("Date of render: ",today)) 
  
  # Select unique top10 genes resolution 0.1
  #test2 <- paste0("FindAllMarkersNegbinom.res",res)
  #top <-  integrated@misc[[test2]]  %>% group_by(cluster) %>% top_n(topN, avg_log2FC)
  #utop<-unique(top$gene)
  #Heatmap of the DGE
  #heatmap.list[[test2]] <- DoHeatmap(integrated, features = top$gene) +    
  #  ggtitle(label = "DGE [Negbinom]",subtitle = paste0("Resolution: ",res)) +
  #  labs(caption= paste0("Date of render: ",today)) 
  
  
  plot_name <- paste0("DGE_heatmap_MAST_res",res,"_",today,".png")
  ggsave(filename= plot_name , plot=heatmap.list[[test1]], path=output.path.plots.DGE, width = 20, height = 20, dpi = 200)
  #plot_name <- paste0("DGE_heatmap_Negbinom_res",res,"_",today,".png")
  #ggsave(filename= plot_name , plot=heatmap.list[[test2]], path=output.path.plots.DGE, width = 20, height = 20, dpi = 200) 
}
########################     DGE ANALYSIS         ##############################


