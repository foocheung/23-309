library(docopt)
library(future)

options(future.globals.maxSize = 64000 * 1024^4)
library("Seurat")
library("dplyr")
library("matrixStats")
library('tidyverse')
library("tidyseurat",lib="./lib") ### GET
library(ggplot2)
library(cowplot)
library(patchwork)
#library(sva)
library(harmony)
library("dsb", lib="./lib") ## GET
library(docopt)
library(future)
library(devtools)
library(SingleCellExperiment)
library(SingleR)
library(celldex)
library(irlba)
library(Matrix)
library(SeuratDisk)

library(tidyverse)
library(symphony, lib="./lib")
library(ggpubr)
library(patchwork)
library(BoneMarrowMap, lib="./lib")

source("dev.23_309_functions_cite.R")

process_dataset <- function(config_file) {
  # Load configuration from YAML file
  config <- yaml::read_yaml(config_file)
  
  # B1_US_data <- Read10X_h5(config$data_path)
  B1_US_data <- read_h5_files(config$lanes, config$rawdir, config$prefix)
  
  B1_US_SeuratObj <- create_seurat_objects(B1_US_data, config$lanes)

  B1_US_SeuratObj <- normalize_rna_data(B1_US_SeuratObj) 
  
  B1_US_SeuratObj <- add_metadata_and_filter(B1_US_SeuratObj,config$output_prefix)
 
#  B1_US_SeuratObj<-merge_seurat_objects(B1_US_SeuratObj,config$lanes)

bm <- NormalizeData(B1_US_SeuratObj, assay= "HTO", normalization.method = "CLR", margin=2)
bm <- ScaleData(bm, assay = "HTO", model.use = "negbinom")
bm <- HTODemux(bm, assay = "HTO",positive.quantile = 0.9999)

t<-table(bm$Lane,bm$hash.ID)
write.table(t,paste0(config$output_prefix, "HTO_Lanes.txt",sep=""))

saveRDS(bm, paste0(config$output_prefix, "bm.rds"))
 bm[["RNA"]] <- as(object = bm[["RNA"]], Class = "Assay")

bm <-runTransferLearningBM(bm)
 ggsave(paste0(config$output_prefix,"bone_marrow","_umap.pdf"), width = 26, height = 9, plot = DimPlot(subset(bm, mapping_error_QC == 'Pass'), reduction = 'umap', group.by = c('predicted_CellType'), split.by = "Lane", raster=FALSE, label=TRUE, label.size = 4) & NoLegend()
) 
  
}

process_dataset("config_2.yaml")



