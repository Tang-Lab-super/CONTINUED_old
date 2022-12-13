# /home/yuchen/miniconda3/envs/R4.0/bin/R
# filter m/z max intensity is less than 400

########################################################################################################################
rm(list=ls())

work_dir <- csv_output_dir_path
setwd(work_dir)

"""
制作一个sample.list.txt文件，储存所有的样品名，每行一个sample
"""
sample_list <- read.table(file="sample.list.txt", header=F)$V1

################################################################
for (sample_name in sample_list ){
print(sample_name)
setwd(paste0(work_dir, "/", sample_name))

#把原始文件名改成unfilter后缀
file.rename(paste0(sample_name,"_3k_signal.lock_mass.txt"), paste0(sample_name,"_3k_signal.lock_mass.unfiltered.txt"))
file.rename(paste0(sample_name,"_3k_signal.unlock_mass.txt"), paste0(sample_name,"_3k_signal.unlock_mass.unfiltered.txt"))

signal_total_lock <- read.table(paste0(sample_name,"_3k_signal.lock_mass.unfiltered.txt"), header=F)
signal_total_unlock <- read.table(paste0(sample_name,"_3k_signal.unlock_mass.unfiltered.txt"), header=F)

signal_total_lock_mat <- signal_total_lock[2:nrow(signal_total_lock), 4:ncol(signal_total_lock)]
signal_total_unlock_mat <- signal_total_unlock[2:nrow(signal_total_unlock), 4:ncol(signal_total_unlock)]

rownames(signal_total_lock_mat) <- signal_total_lock[2:nrow(signal_total_lock), 1]
colnames(signal_total_lock_mat) <- signal_total_lock[1, 4:ncol(signal_total_lock)]

rownames(signal_total_unlock_mat) <- signal_total_unlock[2:nrow(signal_total_unlock), 1]
colnames(signal_total_unlock_mat) <- signal_total_unlock[1, 4:ncol(signal_total_unlock)]

signal_total_lock_mat <- signal_total_lock_mat[,apply(signal_total_lock_mat, 2, max) >= 400]
signal_total_unlock_mat <- signal_total_unlock_mat[,apply(signal_total_unlock_mat, 2, max) >= 400]

signal_total_lock_mat_sel_names <- c("Index", "x", "y", colnames(signal_total_lock_mat))
signal_total_unlock_mat_sel_names <- c("Index", "x", "y", colnames(signal_total_unlock_mat))

signal_total_lock_sel <- signal_total_lock[, signal_total_lock[1,] %in% signal_total_lock_mat_sel_names]
signal_total_unlock_sel <- signal_total_unlock[, signal_total_unlock[1,] %in% signal_total_unlock_mat_sel_names]

write.table(signal_total_lock_sel, file=paste0(sample_name,"_3k_signal.lock_mass.txt"), quote=F, sep="\t", row.names=F, col.names=F)
write.table(signal_total_unlock_sel, file=paste0(sample_name,"_3k_signal.unlock_mass.txt"), quote=F, sep="\t", row.names=F, col.names=F)
}
########################################################################################################################
