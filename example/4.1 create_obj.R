#Create desi object with best npcs and resolution.

sample_name <- 'ST87_20210331'
nor_method <- 'ScaleData'
npcs <- 40
res <- 0.6
setwd(csv_output_dir_path)
########################################################################################################################
# unlock mass
input_unlockmass <- paste0(sample_name, "_signal.unlock_mass.txt")
signal_total_unlock <- read.table(input_unlockmass, header=F)
#signal_total_unlock[2:length(signal_total_unlock$V2),]$V2 <- as.character(as.numeric(signal_total_unlock[2:length(signal_total_unlock$V2),]$V2))
#signal_total_unlock[2:length(signal_total_unlock$V3),]$V3 <- as.character(as.numeric(signal_total_unlock[2:length(signal_total_unlock$V3),]$V3))
mass_selected_unlock <- as.data.frame(as.numeric(signal_total_unlock[1, 4:ncol(signal_total_unlock)]))
colnames(mass_selected_unlock) <- "M.Z"

# lock mass
input_lockmass <- paste0(sample_name, "_signal.lock_mass.txt")
signal_total_lock <- read.table(input_lockmass, header=F)
#signal_total_lock[2:length(signal_total_lock$V2),]$V2 <- as.character(as.numeric(signal_total_lock[2:length(signal_total_lock$V2),]$V2))
#signal_total_lock[2:length(signal_total_lock$V3),]$V3 <- as.character(as.numeric(signal_total_lock[2:length(signal_total_lock$V3),]$V3))
mass_selected_lock <- as.data.frame(as.numeric(signal_total_lock[1, 4:ncol(signal_total_lock)]))
colnames(mass_selected_lock) <- "M.Z"


# coordinates covered by tissue
input_coor <- paste0(sample_name, "_signal.selected.coor.txt")
coor_selected <- read.table(input_coor, header=T)
coor_selected$x <- as.character(coor_selected$x)
coor_selected$y <- as.character(coor_selected$y)
########################################################################################################################

desi_lock_mass_obj <- create_desi_obj_with_res(signal_total_lock, mass_selected_lock, coor_selected, npcs, nor_method, res)
desi_unlock_mass_obj <- create_desi_obj_with_res(signal_total_unlock, mass_selected_unlock, coor_selected, npcs, nor_method, res)

desi_lock_mass_obj_adjust <- desi_lock_mass_obj
desi_lock_mass_obj_adjust@meta.data[, c("seurat_clusters")] <- NA
desi_lock_mass_obj_adjust@meta.data[rownames(desi_lock_mass_obj_adjust@meta.data) %in% rownames(desi_unlock_mass_obj@meta.data), ]$seurat_clusters <- desi_unlock_mass_obj@meta.data$seurat_clusters
Idents(desi_lock_mass_obj_adjust) <- desi_lock_mass_obj_adjust@meta.data$seurat_clusters

###########################################
dir <- paste0(csv_output_dir_path, '/', sample_name, '/', 'run.nor_', nor_method, '.npcs_', npcs, '.res_', res)
dir.create(paste0(csv_output_dir_path, '/', sample_name))
dir.create(dir)

output_rds <- paste0(dir, "/", sample_name, ".desi.lock.mass.rds")
saveRDS(desi_lock_mass_obj, file=output_rds)
rm(output_rds)

output_rds <- paste0(dir, "/", sample_name, ".desi.unlock.mass.rds")
saveRDS(desi_unlock_mass_obj, file=output_rds)
rm(output_rds)

output_rds <- paste0(dir, "/", sample_name, ".desi.lock.mass.adjust.rds")
saveRDS(desi_lock_mass_obj_adjust, file=output_rds)
rm(output_rds)
###########################################
plot_seurat_clusters(desi_lock_mass_obj, signal_total_lock, paste0(npcs, "_", res, "_lockmass"), dir, col)
plot_seurat_clusters(desi_unlock_mass_obj, signal_total_unlock, paste0(npcs, "_", res, "_unlockmass"), dir, col)
coor_border_cluster_all <- plot_seurat_clusters_border(desi_unlock_mass_obj, signal_total_unlock, paste0(npcs, "_", res, "_unlockmass_border"), dir)

write.table(coor_border_cluster_all, file=paste0(dir, "/", "coor_border_cluster_all.txt"), quote=F, sep="\t", row.names=F, col.names=T)
rm(dir)

########################################################################################################################
