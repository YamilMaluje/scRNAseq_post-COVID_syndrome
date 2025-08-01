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
today <- Sys.Date() %>% format(., format="%d.%m.%Y")
output.path <- paste(orig,"05_DownstreamAnalysis","SelectedGenes",today,sep = "/")
output.path.plots <- paste(output.path, "Plots",sep = "/")
dir.create(output.path.plots,recursive = TRUE)
clu <- c(0.3)
proj.name <- "AFL_NEW_All"
specie<- "hsa" # homo sapiens "hsa" or mus musculus "mmu"
################# SETTING VARIABLES ##############################################
#######################    FUNCTIONS      ######################################
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
  x
}
#######################    FUNCTIONS      ######################################
integrated <- readRDS("/Volumes/MarkV/Projects/UKSH scRNAseq/Server Files/scRNAseq Analysis/COVID/COVID-AFL_NEW_ALL/05_DownstreamAnalysis/09.04.2025/AFL_NEW_All.integrated_od_markers_annotated_09.04.2025.rds")


################# ADDING METADATA TO SEURAT INTEGRATED OBJECT  #################
# Create a tibble (dataframe) for the cluster data
cluster_data <- tibble(
  Cluster = c(1, 2, 3, 4, 5, 6),
  Cell_Type = c("Differentiating Basal cells",
                "Secretory Cell",
                "Basal Cells",
                "EMT Basal-like cells",
                "Ciliated Cells" ,
                "EMT Epithelial cells"),
  Markers = list(
    c("KRT23","KRT6B", "SPDEF", "KRT7"),
    c("CYP2F1","SCGB1A1", "SCGB3A1" ),
    c("KRT5", "TP63", "TP63","KRT6A", "KRT13"),
    c("VIM", "KRT6A","FAP","LUM"),
    c("TUBB4B","DNAI1","FOXJ1", "RSPH4A"),
    c("FAP", "COL16A1","VIM","DCN")
  )
)



# Convert the list of vectors in the Markers column into comma-separated strings
cluster_data_excel <- cluster_data %>%
  mutate(Markers = sapply(Markers, function(x) paste(x, collapse = ", ")))
# Save the tibble as an Excel file
writexl::write_xlsx(cluster_data_excel,paste(output.path.plots,"CellType_CanonicalMarkers_Annotation_res0.1.xlsx",sep = "/"))



#ADD METADATA ANNOTATION TO SEURAT OBJECT
iden <-  "integrated.clusters.res.0.1.harmony.Batch"
Idents(integrated) <-iden
integrated@meta.data[["CanonicalMarkers_Annotation"]] <- "NULL"
#Go through all of the clusters
for (clus in cluster_data$Cluster){
  message("Old Cluster: ",clus)
  #Subset the metadata to the respective sample
  filter.meta <- subset(cluster_data, Cluster %in% clus)
  #Go trhough the entire metadata column to add it to the seurat object
  integrated@meta.data[["CanonicalMarkers_Annotation"]][integrated@meta.data[[iden]] ==clus] <- filter.meta$Cell_Type
  
}
saveRDS(integrated, file = paste(output.path, paste0("AFL_NEW_All.integrated_od_markers_annotated_",today,".rds"),sep = "/"))
p <- DimPlot(integrated, reduction = "umap.integrated.harmony.Batch" , group.by = "CanonicalMarkers_Annotation", label = T)
save_plot(path =output.path.plots, 
          plot_name = paste0("Umap_CanonicalMarkers_Annotation_",today), plot =p)
################# ADDING METADATA TO SEURAT INTEGRATED OBJECT  #################
