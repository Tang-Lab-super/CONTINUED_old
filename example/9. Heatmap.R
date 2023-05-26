# /home/yuchen/miniconda3/envs/R4.0/bin/R

library("harmony")
library("cowplot")
library("pheatmap")
library("Seurat")
library("RColorBrewer")
library("ggplot2")
library("glue")

########################################################################################################################
infor <- as.data.frame(read.table(file="/data/bingling/rPackageTutorial/Continued/example/Dataset/file.infor.txt", header=T))
dir_path <- '/data/tang/desi/processing/supplementary/'   #obj_path的上级目录
dir_mass_merged <- '/data/tang/desi/processing/merge.lockmass/merged.mass.for.sample.cytof/'   #'{sample_name}_mass.txt'的上级目录

for (l in 1:dim(infor)[1]){
	sample_name <- infor[l, 1]
	obj_path <- infor[l, 2]
	assign(paste0('obj_', sample_name), Create_New_Obj(sample_name, dir_path, obj_path, dir_mass_merged))

	annotated_mass <- read.table(file="/data/bingling/rPackageTutorial/Continued/example/Dataset/annotated.mass.index.txt", header=T)
	annotated_mass <- paste0(rep('Mass-Index-', length(annotated_mass$Index)), annotated_mass$Index)
	annotated_mass <- sub('_', '-', annotated_mass)

	cluster_list <- as.numeric(unique(Idents(get(paste0('obj_', sample_name)))))

	nrow <- dim(get(paste0('obj_', sample_name))[['RNA']]@data)[1]
	ncol <- length(cluster_list)
	mat <- matrix(nrow=nrow, ncol=ncol)

	for (i in 1:length(cluster_list)){
		ident <- cluster_list[i]
		expr_average <- AverageExpression(subset(get(paste0('obj_', sample_name)), idents=ident), slot='scale.data')
		mat[,i] <- as.vector(expr_average[['RNA']][,'all'])
	}

	rownames(mat) <- rownames(get(paste0('obj_', sample_name))[['RNA']]@scale.data)
	colnames(mat) <- paste0(rep('C_', length(cluster_list)), cluster_list)

	mat_anno_mass <- mat[rownames(mat) %in% annotated_mass, ]
	cols <- rev(brewer.pal(11, "Spectral"))
	pdf(file=paste0(sample_name, ".Heatmap.Annotated_mass.by.clusters.pdf"), width=10, height=10)
	p1 <- pheatmap(mat_anno_mass, cluster_rows=T, cluster_cols=T, clustering_method="ward.D2", color=cols)
	print(p1)
	dev.off()

	pdf(file=paste0(sample_name, ".Heatmap.Annotated_mass.by.clusters.scaled.pdf"), width=10, height=10)
	p3 <- pheatmap(mat_anno_mass, cluster_rows=T, cluster_cols=T, clustering_method="ward.D2", color=cols, scale="row")
	print(p3)
	dev.off()
  #mat_anno_mass[p3$tree_row$order,]把原始排序变成聚类后的排序

	pdf(file=paste0(sample_name, ".Heatmap.ALL.by.clusters.pdf"), width=10, height=10)
	p2 <- pheatmap(mat, cluster_rows=T, cluster_cols=T, clustering_method="ward.D2", color=cols)
	print(p2)
	dev.off()

	pdf(file=paste0(sample_name, ".Heatmap.ALL.by.clusters.scaled.pdf"), width=10, height=10)
	p4 <- pheatmap(mat, cluster_rows=T, cluster_cols=T, clustering_method="ward.D2", color=cols, scale="row")
	print(p4)
	dev.off()

}





