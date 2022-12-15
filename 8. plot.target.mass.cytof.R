# /home/yuchen/miniconda3/envs/R4.0/bin/R

library("harmony")
library("cowplot")
library("pheatmap")
library("Seurat")
library("RColorBrewer")
library("ggplot2")
library("glue")
library(future)
# it is tricky to set options
plan("multiprocess", workers = 30)
options(future.globals.maxSize = 30 * 1024^3)
options(future.seed = TRUE)
options(future.rng.onMisuse="ignore")

########################################################################################################################
rm(list=ls())
setwd(plot.target.mass.path)
source("Load.desi.functions.R")
########################################################################################################################
get_matrix_in_order <- function(df, value_list){
    l <- lapply(value_list, function(x){df[rownames(df)==x, ]})
    return (do.call(rbind,l))
}

Create_New_Obj <- function(sample_name, dir_path, obj_path, dir_mass_merged){
	obj <- readRDS(file=paste0(dir_path, obj_path))

	mass_matrix <- obj[['RNA']]@counts

	mass_merged_tab <- read.table(file=paste0(dir_mass_merged, '/', sample_name, '_mass.txt'), header=T, sep="\t")
	colnames(mass_merged_tab) <- c('Mass_Index', 'Mass_Value')
    #bingling add
    mass_merged_tab$Mass_Value = gsub("0+$","",mass_merged_tab$Mass_Value)
    mass_merged_tab$Mass_Value = gsub("0+,",",",mass_merged_tab$Mass_Value)
  
	mass_merged_tab$Mass_Index <- paste0(rep('Mass-Index-', length(mass_merged_tab$Mass_Index)), mass_merged_tab$Mass_Index)
	mass_merged_tab$Mass_Index <- sub('_', '-', mass_merged_tab$Mass_Index)

	# not good solution for this problem, but currently make the script works
	mass_merged_tab$Mass_Value <- sub('(,[0-9.]+)+', '', mass_merged_tab$Mass_Value)

	mass_merged_tab$Mass_Id <- paste0('mass-', mass_merged_tab$Mass_Value)

	mass_matrix_new <- get_matrix_in_order(mass_matrix, mass_merged_tab$Mass_Id)
	rownames(mass_matrix_new) <- mass_merged_tab$Mass_Index

	meta_anno <- obj@meta.data[,c('index', 'x', 'y')]
	meta_anno$sample_tag <- sample_name

	desi_obj <- CreateSeuratObject(counts=mass_matrix_new, meta.data=meta_anno)
	desi_obj <- RenameCells(desi_obj, add.cell.id=sample_name)
	desi_obj <- ScaleData(desi_obj)
	Idents(desi_obj) <- Idents(obj)
	desi_obj@meta.data$seurat_clusters <- Idents(desi_obj)

	return(desi_obj)
}


########################################################################################################################
"""
制作一个file.infor.txt文件，列名"Sample","Object",每行分别对应每个样品的名称、desi.lock.mass.adjust.rds文件路径
"""

infor <- as.data.frame(read.table(file="file.infor.txt", header=T))
dir_path <- output_dir   #所有样品用最佳参数做聚类生成的rds文件所在路径
"""
从merge.mass.pl的输出文件colon_cancer_desi.clustered_mass.table.with.anno.txt中获取每个样品的sample_mass.txt，
列名"Index","sample_mass"，分别是样品m/z对应索引和m/z，
所有sample_mass.txt存储在merged.mass.for.sample/路径下
"""
dir_mass_merged <- 'merged.mass.for.sample/'

for (l in 1:dim(infor)[1]){
	sample_name <- infor[l, 1]
	obj_path <- infor[l, 2]
	assign(paste0('obj_', sample_name), Create_New_Obj(sample_name, dir_path, obj_path, dir_mass_merged))
"""
从merge.mass.pl的输出文件colon_cancer_desi.clustered_mass.table.with.anno.txt中获取第一列，作为annotated.mass.index.txt
"""
	annotated_mass <- read.table(file="annotated.mass.index.txt", header=T)
	annotated_mass <- paste0(rep('Mass-Index-', length(annotated_mass$Index)), annotated_mass$Index)
	annotated_mass <- sub('_', '-', annotated_mass)

	intersect_mass <- intersect(annotated_mass, rownames(get(paste0('obj_', sample_name))))

	input_signal_total <- paste0(csv_output_dir_path, '/', sample_name, '/', sample_name, "_3k_signal.lock_mass.txt")
	signal_total <- read.table(input_signal_total, header=F)

	dir_name <- paste0('./', sample_name)
	dir.create(dir_name)
	mass_tag <- 'Mass_clustered'

	for (i in 1:length(intersect_mass)){
		mass_value <- intersect_mass[i]
		plot_MASS_from_desi_obj_RNA(get(paste0('obj_', sample_name)), signal_total, mass_value, mass_tag, dir_name)

	}
}

########################################################################################################################




