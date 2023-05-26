# /home/yuchen/miniconda3/envs/R4.0/bin/R


# it is tricky to set options
plan("multiprocess", workers = 30)
options(future.globals.maxSize = 30 * 1024^3)
options(future.seed = TRUE)
options(future.rng.onMisuse="ignore")

########################################################################################################################
rm(list=ls())
setwd(plot.target.mass.path)
########################################################################################################################
#Create a new file named 'file.infor.txt' with column names of 'Sample' and 'Object'. Each line includes name of each sample,
#'desi.Lock.Mass.Adjust.RDS' file path.
infor <- as.data.frame(read.table(file="/data/bingling/rPackageTutorial/Continued/example/Dataset/file.infor.txt", header=T))
dir_path <- output_dir
#Extract 'Index' and '{sample_name}_mass' for each sample from the 'colon_cancer_desi.clustered_mass.table.with.anno.txt'
#and name as '{sample_name}_mass.txt'. Each sample corresponds to a '{sample_name}_mass.txt', which is
#stored in the dir_mass_merged path.
dir_mass_merged <- 'merged.mass.for.sample/'

for (l in 1:dim(infor)[1]){
	sample_name <- infor[l, 1]
	obj_path <- infor[l, 2]
	assign(paste0('obj_', sample_name), Create_New_Obj(sample_name, dir_path, obj_path, dir_mass_merged))
	#Create a new file named 'annotated.mass.index.txt' with 'Index' from the 'colon_cancer_desi.clustered_mass.table.with.anno.txt'
	annotated_mass <- read.table(file="/data/bingling/rPackageTutorial/Continued/example/Dataset/annotated.mass.index.txt", header=T)
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





