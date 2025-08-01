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
library(ggpubr)
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
integrated <- readRDS("/Volumes/MarkV/Projects/UKSH scRNAseq/Server Files/scRNAseq Analysis/COVID/COVID-AFL_NEW_ALL/05_DownstreamAnalysis/SelectedGenes/09.05.2025/AFL_NEW_All.integrated_od_markers_annotated_09.05.2025.rds")
df <- read_excel_allsheets("/Volumes/MarkV/Projects/UKSH scRNAseq/Server Files/scRNAseq Analysis/COVID/COVID-AFL_NEW_ALL/05_DownstreamAnalysis/09.04.2025/Plots/01_General/FindAllMarkersMAST.res0.1.xlsx")
#Idents(integrated)<- "integrated.clusters.harmony.Sample.res0.3"
library(patchwork)
library(ggpmisc)
library("ggpubr")
iden <- "CanonicalMarkers_Annotation"
redu <- "umap.integrated.harmony.Batch"
#Cluster Annotation markers
canonical_markers <- list("Differentiating Basal cells"=c("KRT23","KRT6B", "SPDEF", "KRT7"),
                          "Secretory Cell"= c("CYP2F1","SCGB1A1", "SCGB3A1" ),
                          "Basal Cells"=c("KRT5", "TP63", "TP63","KRT6A", "KRT13"),
                          "EMT Basal-like cells"=c("VIM", "KRT6A","FAP","LUM"),
                          "Ciliated Cells" =  c("TUBB4B","DNAI1","FOXJ1", "RSPH4A"),
                          "EMT Epithelial cells" =c("FAP", "COL16A1","VIM","DCN"))
################# ThIS IS WORKING BUT NOT PERFECT VERTICAL  ##############################
plot_list <- list() # Lista para almacenar los arreglos de gráficos
# Select unique top10 genes resolution 0.1
for (nm in names(canonical_markers)){
  #Extract the DGE values  of the specific genes
  DGE_val <- df$Sheet1[df$Sheet1$gene %in% canonical_markers[[nm]],] %>% select(., !p_val )
  DGE_val <- DGE_val %>%
    mutate(avg_log2FC = round(avg_log2FC, 2),
           p_val_adj = signif(p_val_adj, 3))  # Adjust to your column names
  
  aux.name <- do.call(paste, c(as.list(canonical_markers[[nm]]), sep = " / "))
  vln_plot <- VlnPlot(integrated,
                      assay = "RNA",
                      features = canonical_markers[[nm]],
                      pt.size = 0.1,
                      raster = FALSE, group.by = iden,stack = T,flip = T) +
    labs(x=NULL,y="Expression (Logcounts)") + NoLegend() +
    theme(aspect.ratio = 0.6, # Ajuste adicional de la relación de aspecto
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          axis.title.y = element_text(size = 9),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) # Reducción adicional de márgenes
  
  
  #This will plot only the combined plot
  DefaultAssay(integrated) <- "RNA"
  Idents(integrated) <- iden
  
  
  # Extract expression values for defining a common scale across all genes for this cell type
  expr_vals <- FetchData(integrated, vars = canonical_markers[[nm]])
  max_expr_all <- max(expr_vals)
  mid_expr_all <- max_expr_all / 2
  
  feature_plots <- lapply(canonical_markers[[nm]], function(gene) {
    FeaturePlot(integrated,
                features = gene,
                raster = FALSE,
                reduction = redu,                              
                label = TRUE,
                repel = TRUE) +
      scale_color_gradientn(
        colors = c("#f0f0f0", "#f2f224", "#d9596a", "#3b0793"),
        values = scales::rescale(c(0, 0.01, mid_expr_all, max_expr_all)),
        limits = c(0, max_expr_all),
        oob = scales::squish
      ) +    theme_void()+
      theme(#aspect.ratio = 1, # Ajustar la relación de aspecto para cada subplot
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #plot.title = element_text(hjust = 0.5, size = 10),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal"
      )# Inicialmente la leyenda está a la derecha
  })
  
  # Extraer la leyenda del primer plot
  #legend <- get_legend(feature_plots[[1]])
  legend = cowplot::get_plot_component(feature_plots[[1]], 'guide-box-bottom', return_all = TRUE)
  
  # Eliminar las leyendas de los plots individuales
  feature_plots_without_legend <- lapply(feature_plots, function(p) p + theme(legend.position = "none"))
  
  # Combinar los plots sin leyenda
  combined_feature_plot <- cowplot::plot_grid(plotlist = feature_plots_without_legend, ncol = 1) # Ajusta ncol según necesites
  
  # Añadir la leyenda común debajo de los plots combinados
  final_feature_plot <- cowplot::plot_grid(combined_feature_plot, legend, ncol = 1, rel_heights = c(length(canonical_markers[[nm]]), 0.3)) & 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  table_plot <- NULL
  if (dim(DGE_val)[1] > 1){
    table_plot <- ggtexttable(DGE_val, rows = NULL,
                              theme = ttheme("blank", base_size = 7)) %>% # Reducción adicional del tamaño del texto base
      tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1) %>% # Líneas más delgadas
      tab_add_hline(at.row = dim(DGE_val)[1]+1, row.side = "bottom", linewidth = 2, linetype = 1) # Línea inferior más delgada
  } else {
    # Crear un gráfico vacío si no hay tabla
    table_plot <- ggplot() + theme_void()
  }
  
  # Combinar los gráficos en una fila de 3 columnas
  combined_plot <- ggarrange(table_plot, final_feature_plot, vln_plot, ncol = 3, widths = c(0.8, 1, 1)) # Ajuste de anchos
  
  # Añadir un título a la página
  combined_plot <- annotate_figure(combined_plot, top = text_grob(paste("Cell Type:", nm), face = "bold", size = 10)) # Reducción del tamaño del título
  
  plot_list[[nm]] <- combined_plot
  
}
plot_name <- paste0("canonical_markers_plots_",today,".pdf")
pdf(file= paste(output.path.plots,plot_name,sep="/"), width = 12, height = 8.27,onefile = TRUE)
for(i in names(plot_list)){
  print(plot_list[[i]])
}
dev.off()
################# ThIS IS WORKING BUT NOT PERFECT VERTICAL ##############################

