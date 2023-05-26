



#' Get desi objects of different resolutions or normalization methods.
#'
#' @param signal_total dataFrame from "_signal.unlock_mass.txt" file
#' @param mass_selected m/z you want
#' @param coor_selected dataFrame from "_signal.selected.coor.txt" file
#' @param npcs number of principal components in RunPCA and FindNeighbors
#' @param nor Normalization methods include SCTransform and ScaleData.
#'
#' @return desi objects of class Seurat
#' @export
#'
#' @examples " "
#'
create_desi_obj <- function(signal_total, mass_selected, coor_selected, npcs, nor){
	signal_total_selected_col_mass <- which(signal_total[1,] %in% mass_selected$M.Z)
	signal_total_coor <- signal_total[2:length(signal_total$V1),1:3]

	colnames(signal_total_coor) <- c("index_total", "x_total", "y_total")
	signal_total_coor <- merge(signal_total_coor, coor_selected, by.x="x_total", by.y="x")
	signal_total_coor_selected <- signal_total_coor[signal_total_coor$y_total == signal_total_coor$y, ]
	meta_coor <- signal_total_coor_selected[,c("index_total", "x_total", "y_total")]
	dim(signal_total_coor_selected)
	dim(coor_selected)
	dim(meta_coor)

	signal_total_selected <- signal_total[signal_total$V1 %in% signal_total_coor_selected$index_total, ]
	meta_cell <- paste0("cell_", signal_total_selected$V1)
	signal_total_selected <- signal_total_selected[,signal_total_selected_col_mass]
	rownames(signal_total_selected) <- meta_cell
	meta_mass <- signal_total[1, signal_total_selected_col_mass]
	colnames(signal_total_selected) <- paste0("mass-", meta_mass)

	meta_anno <- as.data.frame(meta_coor)
	colnames(meta_anno) <- c("index", "x", "y")

	signal_total_selected = as.data.frame(lapply(signal_total_selected,as.numeric))  #bingling add
	rownames(signal_total_selected) <- meta_cell   #bingling add

	meta_anno$intensity_max <- apply(signal_total_selected, 1, max)
	meta_anno$intensity_mean <- apply(signal_total_selected, 1, mean)
	meta_anno$intensity_sum <- apply(signal_total_selected, 1, sum)
	rownames(meta_anno) <- paste0("cell_", meta_anno$index)

	desi_obj <- CreateSeuratObject(counts=t(signal_total_selected), meta.data=meta_anno)

	if(nor=="SCT"){
	desi_obj <- SCTransform(desi_obj)
	} else {
	desi_obj <- ScaleData(desi_obj)
	}

	desi_obj <- FindVariableFeatures(desi_obj)
	desi_obj <- RunPCA(desi_obj, npcs = npcs)
	desi_obj <- FindNeighbors(desi_obj, dims = 1:npcs)
	desi_obj <- FindClusters(desi_obj)
	desi_obj <- FindClusters(desi_obj, resolution = 0.2)
	desi_obj <- FindClusters(desi_obj, resolution = 0.4)
	desi_obj <- FindClusters(desi_obj, resolution = 0.6)
	desi_obj <- FindClusters(desi_obj, resolution = 0.8)
	desi_obj <- FindClusters(desi_obj, resolution = 1.0)
	desi_obj <- FindClusters(desi_obj, resolution = 1.2)
	desi_obj <- FindClusters(desi_obj, resolution = 1.4)
	desi_obj <- FindClusters(desi_obj, resolution = 1.6)
	desi_obj <- FindClusters(desi_obj, resolution = 1.8)
	desi_obj <- RunUMAP(desi_obj, dims = 1:npcs, check_duplicates = FALSE)

#	change cluster index
#	table(desi_obj@meta.data$seurat_clusters)
	desi_obj@meta.data$seurat_clusters <- as.factor(as.numeric(desi_obj@meta.data$seurat_clusters) )
	table(desi_obj@meta.data$seurat_clusters)
	Idents(desi_obj) <- desi_obj@meta.data$seurat_clusters

	return(desi_obj)
}


#' Plot seurat clusters with different resolutions. The normalization method is ScaleData.
#'
#' @param desi_obj the output of 'create_desi_obj'
#' @param signal_total dataFrame from "_signal.unlock_mass.txt" file
#' @param npcs number of principal components in RunPCA and FindNeighbors
#' @param dir_name Set a dir path you want to store the results.
#' @param col type of colors
#'
#' @return Cluster results and dimplots. For each fixed resolution and number
#' of principal components, you get a '.cluster.pdf' and a 'Dimplot.*.pdf'.
#' @export
#'
#' @examples " "
#'
plot_seurat_clusters_with_resolution_RNA <- function(desi_obj, signal_total, npcs, dir_name, col){
  # set clustering resolution: from 0.2 to 1.8 by 0.2 as step
  for (i in seq(0.2, 1.8, 0.2)){
    print(paste0("resolution ", i))

    desi_obj@meta.data$seurat_clusters <- as.factor(as.numeric(desi_obj@meta.data[, paste0("RNA_snn_res.", i)]) )
    Idents(desi_obj) <- desi_obj@meta.data$seurat_clusters

    desi_obj_tissue_df <-data.frame(X=as.numeric(desi_obj@meta.data[, "x"]),
                                    Y=as.numeric(desi_obj@meta.data[, "y"]),
                                    Cluster=as.numeric(desi_obj@meta.data[, "seurat_clusters"]) )

    signal_plot <- signal_total[2:length(signal_total$V2), c(2,3)]
    pixal_x <- length(unique(signal_plot$V2))
    pixal_y <- length(unique(signal_plot$V3))

    plot_x_min <- min(as.numeric(signal_plot$V2))
    plot_x_max <- max(as.numeric(signal_plot$V2))
    plot_y_min <- min(as.numeric(signal_plot$V3))
    plot_y_max <- max(as.numeric(signal_plot$V3))

    plot_frame_df <- data.frame(X=c(plot_x_min, plot_x_min, plot_x_max, plot_x_max),
                                Y=c(plot_y_min, plot_y_max, plot_y_min, plot_y_max),
                                Cluster=rep(max(desi_obj_tissue_df$Cluster)+1, 4) )

    desi_obj_df <- rbind(desi_obj_tissue_df, plot_frame_df)

    col_plot <- col[min(desi_obj_df$Cluster):(max(desi_obj_df$Cluster)-1)]
    col_plot <- c(col_plot, "#f0f0f0")

    print(paste0("npcs: ", npcs, " resolution: ", i))
    output_pdf <- paste0(dir_name, "/", "Tissue.npcs", npcs, ".res_", i, ".cluster.pdf")
    pdf(output_pdf, width=10 * pixal_x / pixal_y, height=10 * pixal_x / pixal_y)
    p2 <- ggplot(data=desi_obj_df, aes(x=X, y=Y, fill=factor(Cluster, levels=seq(min(Cluster):max(Cluster)) ) ) ) +
      geom_tile(aes(width=0.05, height=0.05)) +
      # geom_tile() +
      scale_fill_manual(values = col_plot) +
      theme_bw() +
      theme_void() +
      theme(legend.position = "none") + theme(panel.background=element_blank())
    print(p2)
    dev.off()

    pdf(paste0(dir_name, "/", "DimPlot.cluster.npcs", npcs, ".res_", i, ".pdf"), height=10, width=10)
    p1 <- DimPlot(desi_obj, reduction = "umap", cols=col, label.size = 6, label=TRUE, pt.size=1.2) + NoLegend()
    print(p1)
    dev.off()
  }
}


#' Plot seurat clusters with different resolutions. The normalization method is SCTransform.
#'
#' @param desi_obj the output of 'create_desi_obj'
#' @param signal_total dataFrame from "_signal.unlock_mass.txt" file
#' @param npcs number of principal components in RunPCA and FindNeighbors
#' @param dir_name Set a dir path you want to store the results.
#' @param col type of colors
#'
#' @return Cluster results and dimplots. For each fixed resolution and number
#' of principal components, you get a '.cluster.pdf' and a 'Dimplot.*.pdf'.
#' @export
#'
#' @examples " "
#'
plot_seurat_clusters_with_resolution_SCT <- function(desi_obj, signal_total, npcs, dir_name, col){
  # set clustering resolution: from 0.2 to 1.8 by 0.2 as step
  for (i in seq(0.2, 1.8, 0.2)){
    print(paste0("resolution ", i))

    desi_obj@meta.data$seurat_clusters <- as.factor(as.numeric(desi_obj@meta.data[, paste0("SCT_snn_res.", i)]) )
    Idents(desi_obj) <- desi_obj@meta.data$seurat_clusters

    desi_obj_tissue_df <-data.frame(X=as.numeric(desi_obj@meta.data[, "x"]),
                                    Y=as.numeric(desi_obj@meta.data[, "y"]),
                                    Cluster=as.numeric(desi_obj@meta.data[, "seurat_clusters"]) )

    signal_plot <- signal_total[2:length(signal_total$V2), c(2,3)]
    pixal_x <- length(unique(signal_plot$V2))
    pixal_y <- length(unique(signal_plot$V3))

    plot_x_min <- min(as.numeric(signal_plot$V2))
    plot_x_max <- max(as.numeric(signal_plot$V2))
    plot_y_min <- min(as.numeric(signal_plot$V3))
    plot_y_max <- max(as.numeric(signal_plot$V3))

    plot_frame_df <- data.frame(X=c(plot_x_min, plot_x_min, plot_x_max, plot_x_max),
                                Y=c(plot_y_min, plot_y_max, plot_y_min, plot_y_max),
                                Cluster=rep(max(desi_obj_tissue_df$Cluster)+1, 4) )

    desi_obj_df <- rbind(desi_obj_tissue_df, plot_frame_df)

    col_plot <- col[min(desi_obj_df$Cluster):(max(desi_obj_df$Cluster)-1)]
    col_plot <- c(col_plot, "#f0f0f0")

    print(paste0("npcs: ", npcs, " resolution: ", i))
    output_pdf <- paste0(dir_name, "/", "Tissue.npcs", npcs, ".res_", i, ".cluster.pdf")
    pdf(output_pdf, width=10 * pixal_x / pixal_y, height=10 * pixal_x / pixal_y)
    p2 <- ggplot(data=desi_obj_df, aes(x=X, y=Y, fill=factor(Cluster, levels=seq(min(Cluster):max(Cluster)) ) ) ) +
      geom_tile(aes(width=0.05, height=0.05)) +
      # geom_tile() +
      scale_fill_manual(values = col_plot) +
      theme_bw() +
      theme_void() +
      theme(legend.position = "none") + theme(panel.background=element_blank())
    print(p2)
    dev.off()

    pdf(paste0(dir_name, "/", "DimPlot.cluster.npcs", npcs, ".res_", i, ".pdf"), height=10, width=10)
    p1 <- DimPlot(desi_obj, reduction = "umap", cols=col, label.size = 6, label=TRUE, pt.size=1.2) + NoLegend()
    print(p1)
    dev.off()
  }
}

#' Get desi objects with the optimal number of principal components, normalization methods and resolutions.
#'
#' @param signal_total dataFrame from "_signal.unlock_mass.txt" file
#' @param mass_selected m/z you want
#' @param coor_selected dataFrame from "_signal.selected.coor.txt" file
#' @param npcs,nor,res the optimal number of principal components, normalization
#' methods and resolutions you chose from plot_seurat_clusters_with_resolution_SCT'
#' or 'plot_seurat_clusters_with_resolution_RNA'.
#'
#' @return desi objects of class Seurat
#' @export
#'
#' @examples " "
#'
create_desi_obj_with_res <- function(signal_total, mass_selected, coor_selected, npcs, nor, res){
	signal_total_selected_col_mass <- which(signal_total[1,] %in% mass_selected$M.Z)
	signal_total_coor <- signal_total[2:length(signal_total$V1),1:3]

	colnames(signal_total_coor) <- c("index_total", "x_total", "y_total")
	signal_total_coor <- merge(signal_total_coor, coor_selected, by.x="x_total", by.y="x")
	signal_total_coor_selected <- signal_total_coor[signal_total_coor$y_total == signal_total_coor$y, ]
	meta_coor <- signal_total_coor_selected[,c("index_total", "x_total", "y_total")]
	dim(signal_total_coor_selected)
	dim(coor_selected)
	dim(meta_coor)

	signal_total_selected <- signal_total[signal_total$V1 %in% signal_total_coor_selected$index_total, ]
	meta_cell <- paste0("cell_", signal_total_selected$V1)
	signal_total_selected <- signal_total_selected[,signal_total_selected_col_mass]
	rownames(signal_total_selected) <- meta_cell
	meta_mass <- signal_total[1, signal_total_selected_col_mass]
	colnames(signal_total_selected) <- paste0("mass-", meta_mass)

	meta_anno <- as.data.frame(meta_coor)
	colnames(meta_anno) <- c("index", "x", "y")
	meta_anno$intensity_max <- apply(signal_total_selected, 1, max)
	meta_anno$intensity_mean <- apply(signal_total_selected, 1, mean)
	meta_anno$intensity_sum <- apply(signal_total_selected, 1, sum)
	rownames(meta_anno) <- paste0("cell_", meta_anno$index)

	desi_obj <- CreateSeuratObject(counts=t(signal_total_selected), meta.data=meta_anno)

	if(nor=="SCT"){
	desi_obj <- SCTransform(desi_obj)
	} else {
	desi_obj <- ScaleData(desi_obj)
	}

	desi_obj <- FindVariableFeatures(desi_obj)
	desi_obj <- RunPCA(desi_obj, npcs = npcs)
	desi_obj <- FindNeighbors(desi_obj, dims = 1:npcs)
	desi_obj <- FindClusters(desi_obj)
	desi_obj <- FindClusters(desi_obj, resolution = res)
	desi_obj <- RunUMAP(desi_obj, dims = 1:npcs, check_duplicates = FALSE)

#	change cluster index
#	table(desi_obj@meta.data$seurat_clusters)
	desi_obj@meta.data$seurat_clusters <- as.factor(as.numeric(desi_obj@meta.data$seurat_clusters) )
	table(desi_obj@meta.data$seurat_clusters)
	Idents(desi_obj) <- desi_obj@meta.data$seurat_clusters

	return(desi_obj)
}


#' Plot reconstructed tissue patterns
#'
#' @param desi_obj the output of 'create_desi_obj_with_res'
#' @param signal_total dataFrame from "_signal.unlock_mass.txt" file
#' @param npcs number of principal components in RunPCA and FindNeighbors
#' @param dir_name Set a dir path you want to store the results.
#' @param col type of colors
#'
#' @return reconstructed tissue patterns of all clusters and each cluster
#' @export
#'
#' @examples " "
#'
plot_seurat_clusters <- function(desi_obj, signal_total, npcs, dir_name, col){
	# dir_name <- "seurat_clusters.sct"
	dir.create(dir_name)

	pdf(paste0(dir_name, "/", "DimPlot.cluster.npcs", npcs, ".pdf"), height=10, width=10)
	p1 <- DimPlot(desi_obj, reduction = "umap", cols=col, label.size = 6, label=TRUE, pt.size=1.2) + NoLegend()
	print(p1)
	dev.off()

	# plot reconstructed tissue patterns
	# pre matrix
	desi_obj_tissue_df <-data.frame(X=as.numeric(desi_obj@meta.data[, "x"]),
							 Y=as.numeric(desi_obj@meta.data[, "y"]),
							 Cluster=as.numeric(desi_obj@meta.data[, "seurat_clusters"]) )

	signal_plot <- signal_total[2:length(signal_total$V2), c(2,3)]
	pixal_x <- length(unique(signal_plot$V2))
	pixal_y <- length(unique(signal_plot$V3))

	plot_x_min <- min(as.numeric(signal_plot$V2))
	plot_x_max <- max(as.numeric(signal_plot$V2))
	plot_y_min <- min(as.numeric(signal_plot$V3))
	plot_y_max <- max(as.numeric(signal_plot$V3))

	plot_frame_df <- data.frame(X=c(plot_x_min, plot_x_min, plot_x_max, plot_x_max),
								Y=c(plot_y_min, plot_y_max, plot_y_min, plot_y_max),
								Cluster=rep(max(desi_obj_tissue_df$Cluster)+1, 4) )

	desi_obj_df <- rbind(desi_obj_tissue_df, plot_frame_df)

	# plot all clusters in one
	col_plot <- col[min(desi_obj_df$Cluster):(max(desi_obj_df$Cluster)-1)]
	col_plot <- c(col_plot, "#f0f0f0")

	output_pdf <- paste0(dir_name, "/", "Tissue.npcs", npcs, ".cluster.tile_style.pdf")
	pdf(output_pdf, width=10 * pixal_x / pixal_y, height=10 * pixal_x / pixal_y)
	p2 <- ggplot(data=desi_obj_df, aes(x=X, y=Y, fill=factor(Cluster, levels=seq(min(Cluster):max(Cluster)) ) ) ) +
			geom_tile(aes(width=0.05, height=0.05)) +
			# geom_tile() +
			scale_fill_manual(values = col_plot) +
			theme_bw() +
			theme_void() +
			theme(legend.position = "none") + theme(panel.background=element_blank())
	print(p2)
	dev.off()

	output_pdf <- paste0(dir_name, "/", "Tissue.npcs", npcs, ".cluster.raster_style.pdf")
	pdf(output_pdf, width=10 * pixal_x / pixal_y, height=10 * pixal_x / pixal_y)
	p3 <- ggplot(data=desi_obj_df, aes(x=X, y=Y, fill=factor(Cluster, levels=seq(min(Cluster):max(Cluster)) ) ) ) +
			# geom_tile(aes(width=0.05, height=0.05)) +
			geom_raster() +
			scale_fill_manual(values = col_plot) +
			theme_bw() +
			theme_void() +
			theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background=element_blank())
	print(p3)
	dev.off()

	# plot each cluster
	for (i in seq(min(desi_obj_df$Cluster):((max(desi_obj_df$Cluster)) - 1)) ){
		print(paste0("plot cluster_", i))

		col_cluster <- rep("#f0f0f0", (max(desi_obj_df$Cluster) - min(desi_obj_df$Cluster) ) + 1 )
		col_cluster[i] <- col[i]

		output_pdf <- paste0(dir_name, "/", "Tissue.npcs", npcs, ".cluster_", i, ".pdf")
		pdf(output_pdf, width=10 * pixal_x / pixal_y, height=10 * pixal_x / pixal_y)

		p4 <- ggplot(data=desi_obj_df, aes(x=X, y=Y, fill=factor(Cluster, levels=seq(min(Cluster):max(Cluster)) ) ) ) +
				geom_tile(aes(width=0.05, height=0.05)) +
				scale_fill_manual(values = col_cluster) +
				theme_bw() +
				theme_void() +
				theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background=element_blank())
		print(p4)
		dev.off()
	}

}


#' Plot the border of reconstructed tissue patterns
#'
#' @param desi_obj the output of 'create_desi_obj_with_res'
#' @param signal_total dataFrame from "_signal.unlock_mass.txt" file
#' @param npcs number of principal components in RunPCA and FindNeighbors
#' @param dir_name Set a dir path you want to store the results.
#'
#' @return border of reconstructed tissue patterns of all clusters and each cluster
#' @export
#'
#' @examples " "
#'
plot_seurat_clusters_border <- function(desi_obj, signal_total, npcs, dir_name){

  #dir.create(dir_name)

  # plot reconstructed tissue patterns
  # pre matrix
  desi_obj_tissue_df <-data.frame(X=as.numeric(desi_obj@meta.data[, "x"]),
                                  Y=as.numeric(desi_obj@meta.data[, "y"]),
                                  Cluster=as.numeric(desi_obj@meta.data[, "seurat_clusters"]) )

  signal_plot <- signal_total[2:length(signal_total$V2), c(2,3)]
  pixal_x <- length(unique(signal_plot$V2))
  pixal_y <- length(unique(signal_plot$V3))

  plot_x_min <- min(as.numeric(signal_plot$V2))
  plot_x_max <- max(as.numeric(signal_plot$V2))
  plot_y_min <- min(as.numeric(signal_plot$V3))
  plot_y_max <- max(as.numeric(signal_plot$V3))

  plot_frame_df <- data.frame(X=c(plot_x_min, plot_x_min, plot_x_max, plot_x_max),
                              Y=c(plot_y_min, plot_y_max, plot_y_min, plot_y_max),
                              Cluster=rep(max(desi_obj_tissue_df$Cluster)+1, 4) )

  desi_obj_df <- rbind(desi_obj_tissue_df, plot_frame_df)

  # plot all clusters in one
  col_plot <- col[min(desi_obj_df$Cluster):(max(desi_obj_df$Cluster)-1)]
  col_plot <- c(col_plot, "#f0f0f0")


  # plot border
  print(paste0("finding borders...."))
  mat_file <-  paste0(dir_name, "/", "Tissue.npcs", npcs, ".coordinate_cluster.txt")
  write.table(desi_obj_df, mat_file, sep='\t', quote=F)
  # border_dir <-
  #command = glue('python3 /data/bingling/rPackageTutorial/Continued/R/get_cluster_border.py -i {mat_file} -o {dir_name} -n 10')
  #system(command)

  coor_border_cluster <- data.frame(X=double(), Y=double(), Cluster=character())
  for (i in min(desi_obj_df$Cluster):(max(desi_obj_df$Cluster)-1)){
    input_mat <- paste0(dir_name, "/Cluster_", i, ".txt")
    mat <- as.data.frame(read.table(file=input_mat, sep=","))

    mat_x <- data.frame(X=as.numeric(mat[2:nrow(mat), 1]))
    mat_x$X_index <- 1:(nrow(mat)-1)
    mat_y <- data.frame(Y=as.numeric(mat[1, 2:ncol(mat)]))
    mat_y$Y_index <- 1:(ncol(mat)-1)

    tab <- mat[2:nrow(mat), 2:ncol(mat)]
    tab_x_y <- which(tab!=0, arr.ind=TRUE)
    colnames(tab_x_y) <- c("X_index", "Y_index")
    tab_x_y <- as.data.frame(tab_x_y)

    tab_x_y <- merge(tab_x_y, mat_x, by.x="X_index", by.y="X_index")
    tab_x_y <- merge(tab_x_y, mat_y, by.x="Y_index", by.y="Y_index")

    coor_border <- data.frame(X=tab_x_y$X, Y=tab_x_y$Y)
    coor_border$Cluster <- i
    coor_border_cluster <- rbind(coor_border_cluster, coor_border)
  }


  coor_border_cluster_all <- rbind(coor_border_cluster, plot_frame_df)

  output_pdf <- paste0(dir_name, "/", "Tissue.npcs", npcs, ".border.ALL.pdf")
  pdf(output_pdf, width=10 * pixal_x / pixal_y, height=10 * pixal_x / pixal_y)
  p2 <- ggplot(data=coor_border_cluster_all, aes(x=X, y=Y, fill=factor(Cluster, levels=seq(min(Cluster):max(Cluster)) ) ) ) +
    geom_tile(aes(width=0.05, height=0.05)) +
    scale_fill_manual(values = col_plot) +
    theme_bw() +
    theme_void() +
    theme(legend.position = "none", panel.background=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p2)
  dev.off()

  output_tab <- paste0(dir_name, "/", "Tissue.npcs", npcs, ".border.ALL.txt")
  write.table(coor_border_cluster_all, file=output_tab, quote=F, sep="\t", row.names=F)

  for (i in min(coor_border_cluster_all$Cluster):max(coor_border_cluster_all$Cluster)){
    print(paste0("plot cluster border: ", i))

    coor_border_cluster_sub <- subset(coor_border_cluster, Cluster==i)
    coor_border_cluster_sub <- rbind(coor_border_cluster_sub, plot_frame_df)
    col_plot_sub <- c(col[i], "#f0f0f0")

    output_pdf <- paste0(dir_name, "/", "Tissue.npcs", npcs, ".border.Cluster_", i ,".pdf")
    pdf(output_pdf, width=10 * pixal_x / pixal_y, height=10 * pixal_x / pixal_y)
    p2 <- ggplot(data=coor_border_cluster_sub, aes(x=X, y=Y, fill=factor(Cluster, levels=unique(Cluster) ) ) ) +
      geom_tile(aes(width=0.05, height=0.05)) +
      scale_fill_manual(values = col_plot_sub) +
      theme_bw() +
      theme_void() +
      theme(legend.position = "none", panel.background=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(p2)
    dev.off()
  }
  unlink(paste0(dir_name, "/Cluster_*.png"))
  unlink(paste0(dir_name, "/Cluster_*.txt"))

  return(coor_border_cluster_all)
}



#' Plot the distribution of m/z that you are interested in.
#'
#' @param desi_obj the output of 'create_desi_obj_with_res'
#' @param signal_total dataFrame from "_signal.unlock_mass.txt" file
#' @param mass_value m/z value that you are interested in
#' @param mass_tag Set a label that attached to the file name
#' @param dir_name Set a dir path you want to store the results.
#'
#' @return distribution plot of m/z that you are interested in
#' @export
#'
#' @examples " "
#'

plot_MASS_from_desi_obj_RNA <- function(desi_obj, signal_total, mass_value, mass_tag, dir_name){
	print(paste0("processing mass: ", mass_value))

	desi_obj_df <-data.frame(X=as.numeric(desi_obj@meta.data[, "x"]),
						 	 Y=as.numeric(desi_obj@meta.data[, "y"]),
						 	 Cell_ID=rownames(desi_obj@meta.data),
						 	 Cluster=as.numeric(desi_obj@meta.data[, "seurat_clusters"]) )

	mass_value_tab <- desi_obj[['RNA']]@counts[rownames(desi_obj[['RNA']]@counts) %in% mass_value, ]
	mass_value_df <- data.frame(Intensity=as.vector(mass_value_tab), Cell_ID=names(mass_value_tab))
	rm(mass_value_tab)
	desi_obj_df$Intensity[desi_obj_df$Cell_ID %in% mass_value_df$Cell_ID] <- mass_value_df$Intensity

	intensity_quantile <- quantile(desi_obj_df$Intensity, probs = seq(0, 1, 0.05))
	# below 10%
	desi_obj_df$Intensity[desi_obj_df$Intensity < intensity_quantile[[1]]] <- intensity_quantile[[1]]
	# up 95%
	desi_obj_df$Intensity[desi_obj_df$Intensity > intensity_quantile[[20]]] <- intensity_quantile[[20]]


	signal_total_df <- data.frame(X=signal_total[2:length(signal_total$V2), ]$V2,
								  Y=signal_total[2:length(signal_total$V3), ]$V3,
								  Cell_ID=paste0("cell_", signal_total[2:length(signal_total$V1), ]$V1),
								  Cluster=rep((max(desi_obj_df$Cluster)+1), (length(signal_total$V1)-1)) )
	# signal_total_df$Intensity <- 0
	signal_total_df$Intensity <- NA
	signal_total_df <- signal_total_df[!signal_total_df$Cell_ID %in% desi_obj_df$Cell_ID, ]

	signal_plot <- signal_total[2:length(signal_total$V2), c(2,3)]
	pixal_x <- length(unique(signal_plot$V2))
	pixal_y <- length(unique(signal_plot$V3))
	plot_x_min <- min(as.numeric(signal_plot$V2))
	plot_x_max <- max(as.numeric(signal_plot$V2))
	plot_y_min <- min(as.numeric(signal_plot$V3))
	plot_y_max <- max(as.numeric(signal_plot$V3))

	plot_frame_df <- data.frame(X=c(plot_x_min, plot_x_min, plot_x_max, plot_x_max),
								Y=c(plot_y_min, plot_y_max, plot_y_min, plot_y_max),
								Cell_ID=c("cell_BottomLeft", "cell_TopLeftLeft", "cell_BottomRight", "cell_TopRight"),
								Cluster=rep(max(desi_obj_df$Cluster)+1, 4),
								Intensity=c(NA, NA, NA, NA) )
								# Intensity=c(0, 0, 0, 0) )

	desi_obj_df <- rbind(desi_obj_df,signal_total_df, plot_frame_df)

	p <- ggplot(data=desi_obj_df, aes(x=as.numeric(X), y=as.numeric(Y), fill=Intensity)) +
			geom_raster(hjust = 0, vjust = 0,  interpolate=T) +
			scale_fill_viridis_c(na.value = NA) +
			# scale_colour_viridis_c(na.value = NA) +
			# ggtitle(paste0(plot_title)) +
			theme_bw() +
			theme_void() +
			theme(
			    legend.position = c(.95, .05),
			    legend.justification = c("right", "bottom"),
			    legend.box.just = "right",
			    legend.margin = margin(2, 2, 2, 2),
			    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
			    ) + theme(panel.background=element_blank())
			# theme(legend.position = "none")

	output_pdf <- paste0(dir_name, "/", "Plot.", mass_tag, "_", mass_value, ".pdf")
	# title was added in the plot, which requireded to increase height of the plot
	# pdf(output_pdf, width=10 * pixal_x / pixal_y, height=10 * pixal_x / pixal_y * 1.02)
	pdf(output_pdf, width=10 * pixal_x / pixal_y, height=10 * pixal_x / pixal_y)
	print(p)
	dev.off()

	# add legend plot here
	# boxplot plus color range
}



get_matrix_in_order <- function(df, value_list){
  l <- lapply(value_list, function(x){df[rownames(df)==x, ]})
  return (do.call(rbind,l))
}


#' Get desi object from RDS file
#'
#' @param sample_name the file name of the sample you want to obtain
#' @param dir_path home directory where your '.rds' is stored
#' @param obj_path suffix for 'desi.lock.mass.adjust.rds' file path, subdirectory of 'dir_path'
#' @param dir_mass_merged directory name you want to store the results
#'
#' @return desi object of different sample that you are interested in
#' @export
#'
#' @examples " "
#'
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


#' Get metabolites whose annotated type belongs to lipid.
#'
#' @param colon_sub A subset of 'colon_cancer_desi.clustered_mass.table.with.anno.csv'
#'
#' @return Four matrices with different levels of screening.
#'         lipid_T: Lipid annotations are not empty.
#'         lipid_T1: Lipid annotations are not empty and are unique.
#'         lipid_T2: The lipid structure formula is unique.
#'         lipid_T3: Get the structure formula.
#' @export
#'
#' @examples " "
#'
get_lipid <- function(colon_sub){
  lipid_T  =  colon_sub[(colon_sub$anno_lipid!='None')&(colon_sub$anno_small_mol=='None'),'anno_lipid']
  print(paste0("length(lipid_T):",length(lipid_T)))
  #把有唯一映射的脂质结构式挑出来
  lipid_T1 = c()
  for (i in seq(length(lipid_T))){
    t = strsplit(lipid_T[i],'[|]')[[1]]
    if (length(t) == 1){
      lipid_T1 = c(lipid_T1, strsplit(t,';')[[1]][3])
    }
  }
  print(paste0("length(lipid_T1):",length(lipid_T1)))
  #
  lipid_T2 = c()
  for (i in seq(length(lipid_T1))){
    t = strsplit(lipid_T1[i],',')[[1]]
    if (length(t) == 1){
      lipid_T2 = c(lipid_T2, t)
    }
  }
  lipid_T3 = gsub('[-|+].*','',lipid_T2)
  print(paste0("length(lipid_T3):",length(lipid_T3)))
  out = list(lipid_T, lipid_T1, lipid_T2, lipid_T3)
  return(out)
}


#' Get metabolites whose annotated type belongs to molecule.
#'
#' @param colon_sub A subset of 'colon_cancer_desi.clustered_mass.table.with.anno.csv'
#'
#' @return Four matrices with different levels of screening.
#'         mol_T: Molecule annotations are not empty.
#'         mol_T2: The names of annotations that are not all 'NONE'.
#'         mol_T2uni: The annotations is non-empty, non-NONE, and unique.
#' @export
#'
#' @examples " "
#'
get_mol <- function(colon_sub){
  mol_T  =  colon_sub[(colon_sub$anno_lipid=='None')&(colon_sub$anno_small_mol!='None'),'anno_small_mol']
  print(paste0("length(mol_T):",length(mol_T)))
  mol_T2 = list()
  num = 0
  for (i in seq(length(mol_T))){
    mol_T1 = c()
    t = strsplit(mol_T[i],'[|]')[[1]]
    for (j in seq(length(t))){
      tN = strsplit(t[j],';')[[1]]
      mol_T1 = c(mol_T1, tN[length(tN)])
      mol_T1 = unique(mol_T1)
    }
    if (sum(toupper(mol_T1)!='NONE')>0){   #保留不全是None的代谢物
      num = num +1
      mol_T2[num] = list(mol_T1)
    }
  }
  print(paste0("length(mol_T2):",length(mol_T2)))
  #名称唯一
  mol_T2uni = c()
  for (i in seq(length(mol_T2))){
    if (sum(toupper(mol_T2[[i]])!='NONE')==1){       #只有一个对应名称
      mol_T2uni = c(mol_T2uni, mol_T2[[i]][toupper(mol_T2[[i]])!='NONE'])
    }
  }
  print(paste0("length(mol_T2uni):",length(mol_T2uni)))
  out = list(mol_T, mol_T2, mol_T2uni)
  return(out)
}



