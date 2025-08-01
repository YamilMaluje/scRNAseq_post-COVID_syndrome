#TUTORIAL
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
#https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
################################################################################
#library(edgeR)
#library(RColorBrewer)
library(gage)
#library(EnhancedVolcano)
library(PCAtools)
#library(tximport)
#library(annotables) 
library(dplyr)
#library(stringr)
#library(ggalt)
#library(pheatmap)
#library(viridis)
library(ggplot2)
library("writexl")
library(readxl)
library(tibble) #to use enframe() to convert vector in dataframe
library(pathview)
library(msigdf) # https://github.com/ToledoEM/msigdf
library(msigdbr)
library(gageData)
library(stringr)
library(forcats)
library(lubridate)
library(ggh4x)
################################################################################  
parent_path <- getwd()
today <- Sys.Date() %>% format(., format="%d.%m.%Y")
out_path_plots <- paste(parent_path, "05_DownstreamAnalysis" ,"GSEA_analysis_Condition",today,sep = "/")
dir.create(out_path_plots,recursive = TRUE)
species <- "Homo sapiens"                      
################################################################################
#########################     FUNCTIONS              ###########################
################################################################################


dotplot_pathways_by_condition2 <- function (sig_obj,condition=c("Up","Down"),topN=0,path.remove=NULL,cutoff=NULL){
  if(topN == 0){
    topNN =dim(sig_obj$greater)[1]
  }else{topNN=topN}
  #Split the gage output into two dataframes for each conditions and filter the topN pathways
  resup20 = as.data.frame(sig_obj$greater)[1:topNN,1:5]
  resup20$plog <- -log10(resup20$p.val)
  resup20$type = paste0(condition[1])
  resup20$name <- rownames(resup20)
  resup20$col <- "#378805"
  
  if(topN == 0){
    topNN =dim(sig_obj$less)[1]
  }else{topNN=topN}
  resdwn20 = as.data.frame(sig_obj$less)[1:topNN,1:5]
  resdwn20$plog <- -log10(resdwn20$p.val)
  resdwn20$type = paste0(condition[2])
  resdwn20$name <- rownames(resdwn20)
  resdwn20$col <- "#FF0000"
  #Remove the undesire name from the patways
  goall <- rbind(resup20,resdwn20)
  goall$name <- gsub(path.remove,"",goall$name)
  goall$name <- gsub("_"," ",goall$name)
  #goall <-goall[-grep("^NA|^NA.",goall$name),]
  goall <-goall[!grepl("^NA|^NA.",goall$name),]
  #Dotpplot the different pathways of the two conditions
  
  #Replace PBS with unstimulated
  goall$type <- gsub("PBS", "Unstimulated", goall$type)
  # Asignamos los colores basados en los valores de la columna "type"
  goall$back_color[goall$type == "Unstimulated"] <- "#FAFC97"
  goall$back_color[goall$type == "TNFα"] <- "#77DE79"
  goall$back_color[goall$type == "TGFβ"] <- "#FCB249"
  goall$back_color[goall$type == "TNFα_TGFβ"] <- "#B29CDB"
  
  # 1. Create a named vector: names = unique types, values = their unique back_color
  facet_colors <- goall %>%
    distinct(type, back_color) %>%
    deframe()  # This creates a named vector: names = type, values = back_color

  # Get the levels in the order ggplot2 will facet them
  facet_levels <- levels(factor(goall$type))
  
  # Reorder the fill vector accordingly
  fill_vector <- facet_colors[facet_levels]
  
  p <- ggplot(goall, aes(x = plog, y = fct_reorder(name, plog))) +
    geom_point(aes(size = set.size, color = plog)) +
    theme_bw(base_size = 14) +
    scale_colour_gradient(
      limits = c(0, (-log10(goall$p.val) %>% max %>% ceiling)),
      low = "red"
    ) +
    ylab(NULL) +
    xlab("P-value [-log10]") +
    ggh4x::facet_wrap2(
      ~type,
      strip = ggh4x::strip_themed(
        background_x = lapply(fill_vector, function(color) {
          element_rect(fill = color)
        }),
        text_x = list(element_text(color = "black", face = "bold"))
      )
    )+theme(legend.position = "bottom")
  #Barplot horizontal the different pathways of the two conditions
  pl <- ggplot(goall, aes(y = plog, x = fct_reorder(name, plog), fill= col)) +
    geom_bar(stat="identity",show.legend = FALSE) +
    theme_bw(base_size = 14) +
    coord_flip()+
    xlab(NULL) + ylab("P-value [-log10]")+
    ggtitle(paste0(condition[1]," v/s ", condition[2]),
            subtitle = paste0(path.remove," pathways ", "\nGreen: Upregulated | Red: Downregulated","\nqcutoff=",cutoff))+
    labs(caption= paste0("Date of render: ",today)) +
    #scale_fill_manual(values=c("#378805", "#FF0000"))
    scale_fill_identity()
  
  
  plots <-list("dotplot"=p,"barplot"=pl)
  #return the plot
  return(plots)
}
pathway_calculations_by_condition2 <- function (DGE_object=DGE_object,pathway.df=NULL,patway.list=NULL,condition=c("Up","Down"),qcutoff=NULL,topN=0,path.remove=NULL){
  #This function adds the gene list of the pathways
  #gage function requires a named vector of fold changes, where the names of the values are the Entrez gene IDs
  #change the column name to be merge with the entrenez dataframe
  names(DGE_object)[names(DGE_object)=="genes"] <- "symbol"
  #merging both dataset that match the symbol names
  new_df <- merge(x=DGE_object, y=pathway.df, by = "symbol",all.x=TRUE) %>% distinct(., entrez,.keep_all=TRUE)
  #Eliminate all rows that did not found any matching entrenez 
  new_df <- new_df[!is.na(new_df$entrez), ]
  
  #Create the vector with the logfold changes and entrenez names
  foldchanges <- new_df$avg_log2FC
  names(foldchanges) <- new_df$entrez
  
  df_geneLits <- new_df %>%
    group_by(gs_name) %>%
    summarize(symbol_list = paste(symbol, collapse = ","))
  
  #Get results 
  
  DGE_object.rt <- gage(foldchanges, 
                        gsets = patway.list,
                        same.dir = TRUE) %>% 
    sigGeneSet(.,cutoff = qcutoff,qpval = "q.val", heatmap = FALSE)
  lapply(DGE_object.rt,head)
  
  #dotplot the results
  plot <- dotplot_pathways_by_condition2(DGE_object.rt,condition=condition,topN=topN,path.remove=path.remove, cutoff=qcutoff)
  
  #Save pathways into tables. each up/down/stats will be in a different sheet
  great.df <- as.data.frame(DGE_object.rt$greater) %>% cbind("Pathway"=rownames(DGE_object.rt$greater))
  less.df <- as.data.frame(DGE_object.rt$less) %>% cbind("Pathway"=rownames(DGE_object.rt$less))
  stats.df <- as.data.frame(DGE_object.rt$stats) %>% cbind("Pathway"=rownames(DGE_object.rt$stats))
  
  
  # Merge the dataframes based on gs_name and Pathway
  if(dim(great.df)[1] >0) {
    great.df <- merge(great.df,df_geneLits, by.x = "Pathway", by.y = "gs_name", all.x = TRUE)
    names(great.df)[names(great.df) == "symbol_list"] <- "Genes"}
  if(dim(less.df)[1] >0) {
    less.df <- merge(less.df,df_geneLits, by.x = "Pathway", by.y = "gs_name", all.x = TRUE)
    names(less.df)[names(less.df) == "symbol_list"] <- "Genes"}
  if(dim(stats.df)[1] >0) {
    stats.df <- merge(stats.df,df_geneLits, by.x = "Pathway", by.y = "gs_name", all.x = TRUE)
    names(stats.df)[names(stats.df) == "symbol_list"] <- "Genes"}
  
  
  
  obj.list <-list("plots"=plot,"gage"=list("great"=great.df,"less"=less.df,"stats"=stats.df))
  return(obj.list)
  
}
caculate_pathways <- function(){
  #qcutoff <- c(0.05,0.1,0.5)
  qcutoff <- c(0.05)
  
  for (q in qcutoff){
    message("Q cutooff value=",q)

    
    message("HALLMARK")
    #############################       HALLMARK            ########################
    #############################    EpCAM vs Control        #######################
    res.all <-pathway_calculations_by_condition2(DGE_object=df.table.cond,
                                                 pathway.df=h.hallmark.df,
                                                 patway.list=h.hallmark.df.list,
                                                 condition=cond,
                                                 qcutoff=q,
                                                 topN=30,
                                                 path.remove="HALLMARK_")
    
    if (dim(res.all$plots[["dotplot"]]$data)[1] !=0){
      out_path_plot <- paste(out_path_plots,paste0("Pathway_dotplot_HALLMARK_",cond[1],"vs",cond[2],"_qcut",q,".png"), sep = "/")
      ggsave(plot = res.all$plots[["dotplot"]] ,filename = out_path_plot, width = 14, units = "in",dpi = 600) 
      
      out_path_plot <- paste(out_path_plots,paste0("Pathway_barplot_HALLMARK_",cond[1],"vs",cond[2],"_qcut",q,".pdf"), sep = "/")
      ggsave(plot = res.all$plots[["barplot"]] ,filename = out_path_plot, width = 14, units = "in",dpi = 600) 
      #returb excel
      write_xlsx(res.all[["gage"]],
                 paste(out_path_plots,paste0("Pathway_HALLMARK_",cond[1],"vs",cond[2],"_qcut",q,".xlsx"),sep = "/"))
      
    }else(message("No pathways found..."))
    
    
  }
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
################################################################################
#########################     FUNCTIONS              ###########################
################################################################################
# Gene sets ---------------------------------------------------------------
# MSigDB v7.4
c2.rt.df = msigdbr(species = species,category = "C2") %>%
  filter(grepl("*REACTOME", gs_subcat)) %>%
  dplyr::select(gs_name, symbol = gene_symbol, entrez= entrez_gene) %>%
  group_by(gs_name) 
# get all reactome pathways
c2.rt.df.list = msigdbr(species = species,category = "C2") %>%
  filter(grepl("*REACTOME", gs_subcat)) %>%
  dplyr::select(gs_name, symbol = gene_symbol, entrez= entrez_gene) %>%
  group_by(gs_name)  %>%
  summarize(symbol=list(entrez)) %>%
  deframe()
#get all Hallmark pathways
h.hallmark.df <- msigdbr(species = species,category = "H")  %>%
  filter(grepl("*HALLMARK", gs_name)) %>%
  dplyr::select(gs_name, symbol = gene_symbol, entrez= entrez_gene) %>%
  group_by(gs_name) 
h.hallmark.df.list <- msigdbr(species = species,category = "H")  %>%
  filter(grepl("*HALLMARK", gs_name)) %>%
  dplyr::select(gs_name, symbol = gene_symbol, entrez= entrez_gene) %>%
  group_by(gs_name)  %>%
  summarize(symbol=list(entrez)) %>%
  deframe()
################################################################################
#SINGLE CELL OUTPUT APPROACH
#loading the logFC results of the different conditions
df <- read_excel_allsheets(paste(parent_path,
                                 "05_DownstreamAnalysis",
                                 "DGE_scRNAseq_Condition",
                                 "11.04.2025",
                                 "FindMarkersMAST_Condition_comparison.xlsx",sep = "/"))


#Go through all o the different comparisons
for(comp in names(df)){
  cond <- str_split_1(comp, pattern = "_vs_")
  df.table.cond <- df[[comp]]
  
  out_path_plots <- paste(parent_path, "05_DownstreamAnalysis" ,"GSEA_analysis_Condition",today,sep = "/")
  out_path_plots <- paste(out_path_plots, comp,sep = "/")
  dir.create(out_path_plots,recursive = TRUE)
  
  caculate_pathways()
}

