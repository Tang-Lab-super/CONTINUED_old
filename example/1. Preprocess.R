
library(CONTINUED)
library(reticulate)

library("harmony")
library("cowplot")
library("pheatmap")
library("Seurat")
library("RColorBrewer")
library("ggplot2")
library("glue")
library("future")

source_python('/data/bingling/rPackageTutorial/Continued/R/Load.python.functions.py')
repl_python()

#########################################################################################
#run in python
data_file_path = '/data/bingling/rPackageTutorial/Continued/example/Dataset/20210331_ST87un.txt'
data_file_path_lm = '/data/bingling/rPackageTutorial/Continued/example/Dataset/20210331_ST87.txt'
mz_from_path = '/data/bingling/rPackageTutorial/Continued/example/Dataset/'
csv_output_dir_path = '/data/bingling/rPackageTutorial/Continued/example/Dataset/'

ST87_20210331 = pd.read_csv(data_file_path, sep = '\t')
ST87_20210331_lm = pd.read_csv(data_file_path_lm, sep = '\t')   #lockmass
sample='ST87_20210331'

df_data = find_interest_factor(ST87_20210331,
                               mz_from_path,
                               mz='788.5479',
                               threshold=500)

mzbg, mz_tissue, predata_signal = remove_bl(ST87_20210331, df_data, 90, 0.5)
ST87_20210331_signal, ST87_20210331_lm_signal, coor = get_output_csv_file(ST87_20210331, ST87_20210331_lm, df_data, csv_output_dir_path, sample)
quit

#########################################################################################
#go back to R


