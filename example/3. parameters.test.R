# test npcs and find best resolution

sample_name = 'ST87_20210331'
#setwd(paste0(csv_output_dir_path, sample_name))
getwd()
########################################################################################################################
# read files
input_unlockmass <- paste0(sample_name, "_signal.unlock_mass.txt")
#signal_total_unlock <- read.table(input_unlockmass, header=F)  #Tang
signal_total_unlock <- read.table(input_unlockmass, header=T, check.names=F)  #bingling
signal_total_unlock = rbind(colnames(signal_total_unlock),signal_total_unlock)
colnames(signal_total_unlock) = paste0('V',seq(ncol(signal_total_unlock)))

#Tang
#signal_total_unlock[2:length(signal_total_unlock$V2),]$V2 <- as.character(as.numeric(signal_total_unlock[2:length(signal_total_unlock$V2),]$V2))
#signal_total_unlock[2:length(signal_total_unlock$V3),]$V3 <- as.character(as.numeric(signal_total_unlock[2:length(signal_total_unlock$V3),]$V3))

mass_selected_unlock <- as.data.frame(as.numeric(signal_total_unlock[1, 4:ncol(signal_total_unlock)]))
colnames(mass_selected_unlock) <- "M.Z"

# coordinates covered by tissue
input_coor <- paste0(sample_name, "_signal.selected.coor.txt")
coor_selected <- read.table(input_coor, header=T)
coor_selected$x <- as.character(coor_selected$x)
coor_selected$y <- as.character(coor_selected$y)

########################################################################################################################
# npcs: 20, 25, 30, 35, 40

col <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
         '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
         '#fffac8', '#800000', '#aaffc3', '#808000',
         #'#ffd8b1', '#000075', '#808080', '#ffffff', '#000000',
         rev(brewer.pal(n = 11, name = "Spectral")),
         rev(brewer.pal(n = 11, name = "PiYG")),
         # rev(brewer.pal(n = 12, name = "Set3")),
         rev(brewer.pal(n = 8, name = "Accent")) )

nor_method <- "ScaleData"
dir <- paste0("parameter.npcs.and.res.test.", nor_method)
dir.create(dir)

for (npcs in c(20, 25, 30, 35, 40)){
        print(paste0("npcs: ", npcs))
        desi_unlock_mass_obj <- create_desi_obj(signal_total_unlock, mass_selected_unlock, coor_selected, npcs, nor_method)
        plot_seurat_clusters_with_resolution_RNA(desi_unlock_mass_obj, signal_total_unlock, npcs, dir, col)
}
rm(dir)

#nor_method <- "SCT"
#dir <- paste0("parameter.npcs.and.res.test.", nor_method)
#dir.create(dir)
#for (npcs in c(20, 25, 30, 35, 40)){
#        print(paste0("npcs: ", npcs))
#        desi_unlock_mass_obj <- create_desi_obj(signal_total_unlock, mass_selected_unlock, coor_selected, npcs, nor_method)
#        plot_seurat_clusters_with_resolution_SCT(desi_unlock_mass_obj, signal_total_unlock, npcs, dir, col)
#}
#rm(dir)





