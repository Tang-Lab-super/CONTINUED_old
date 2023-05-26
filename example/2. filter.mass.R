# filter m/z max intensity is less than 400

########################################################################################################################
work_dir <- csv_output_dir_path
setwd(work_dir)

#If you want to batch multiple datasets, you can make a file called 'sample.list.txt' to store all sample
#names, one sample per line. For example:
#sample_list <- read.table(file="sample.list.txt", header=F)$V1

################################################################

sample_name = 'ST87_20210331'
#setwd(paste0(work_dir, "/", sample_name))

#Change the original file name to the 'unfiltered' suffix.
file.rename(paste0(sample_name,"_signal.lock_mass.txt"), paste0(sample_name,"_signal.lock_mass.unfiltered.txt"))
file.rename(paste0(sample_name,"_signal.unlock_mass.txt"), paste0(sample_name,"_signal.unlock_mass.unfiltered.txt"))

signal_total_lock <- read.table(paste0(sample_name,"_signal.lock_mass.unfiltered.txt"), header=F)
signal_total_unlock <- read.table(paste0(sample_name,"_signal.unlock_mass.unfiltered.txt"), header=F)

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

signal_total_lock_sel = cbind(signal_total_lock[,1:3], signal_total_lock_sel)
signal_total_unlock_sel = cbind(signal_total_unlock[,1:3], signal_total_unlock_sel)

write.table(signal_total_lock_sel, file=paste0(sample_name,"_signal.lock_mass.txt"), quote=F, sep="\t", row.names=F, col.names=F)
write.table(signal_total_unlock_sel, file=paste0(sample_name,"_signal.unlock_mass.txt"), quote=F, sep="\t", row.names=F, col.names=F)


