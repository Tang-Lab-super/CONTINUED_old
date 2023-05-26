
#Take the tumor region from different samples as an example to make the input file of WGCNA.
args = commandArgs(trailingOnly=TRUE)


sample_name <- args[1]
cluster_target1 <- args[2] # ScaleData or SCT
cluster_target2 <- args[3] # ScaleData or SCT
cluster_target3 <- args[4] # ScaleData or SCT
cluster_target4 <- args[5] # ScaleData or SCT
cluster_target5 <- args[6] # ScaleData or SCT
cluster_target6 <- args[7] # ScaleData or SCT
cluster_target7 <- args[8] # ScaleData or SCT
cluster_target8 <- args[9] # ScaleData or SCT


# read files
input_rds <- paste0(sample_name, ".desi.lock.mass.adjust.rds")
# For example, '/data/bingling/rPackageTutorial/Continued/example/Dataset/ST87_20210331/run.nor_ScaleData.npcs_40.res_0.6/ST87_20210331.desi.lock.mass.adjust.rds'
rds <- readRDS(input_rds)

############################################################################################################
setwd(paste0("/data/tang/desi/input/combined/", sample_name))  #Home directory of all sample input files.
getwd()
print(sample_name)

# read files
input_lc <- paste0(sample_name, "_3k_signal.lock_mass.txt")
lockmass <- read.csv(input_lc,sep = '\t', header=F)

########################################################################################################################
lockmass[2:length(lockmass$V2),]$V2 <- as.character(as.numeric(lockmass[2:length(lockmass$V2),]$V2))
lockmass[2:length(lockmass$V3),]$V3 <- as.character(as.numeric(lockmass[2:length(lockmass$V3),]$V3))
colnames(lockmass) = lockmass[1,]
lockmass = lockmass[-1,]
lockmass = lockmass[,-1]
rownames(lockmass)=seq(nrow(lockmass))
#step1：标准化
quzz <- function(data){
  q1 = quantile(data,0.01)
  q99 = quantile(data,0.99)
  data[data<q1]=q1
  data[data>q99]=q99
  data = (data-min(data))/(max(data)-min(data))
  return(data)
}
rela_ST88 = apply(lockmass[,3:ncol(lockmass)], 2, quzz)
rela_ST88 = cbind(lockmass[,1:2],rela_ST88)
#step2：抽取肿瘤部位
coor_T = rds@meta.data[rds@meta.data$seurat_clusters %in% c(cluster_target1,cluster_target2,cluster_target3,cluster_target4,cluster_target5,cluster_target6,cluster_target7,cluster_target8),c('x','y')]
library(plyr)
rela_ST88_T = ddply(coor_T,.(x,y),function(x){rela_ST88[rela_ST88$x==x[1,1] & rela_ST88$y==x[1,2],]})
rela_ST88_Tcp = rela_ST88_T
print(dim(rela_ST88_T))
print(dim(rela_ST88))
print(dim(coor_T))
#step3：将m/z更换为inedx
colon = read.csv(file="/data/bingling/rPackageTutorial/Continued/example/Dataset/colon_cancer_desi.clustered_mass.table.with.anno.csv", header=T)
#colon$Index = gsub('_','-',colon$Index)
#colon$new_Index = paste0('Mass-Index-',colon$Index)
index_ST88 = colon[,c('Index',paste0(sample_name,'_mass'))]
print(sum(index_ST88[,2]!='None'))
newcol = index_ST88$Index[match(colnames(rela_ST88_T)[3:ncol(rela_ST88_T)],index_ST88[,2])]
colnames(rela_ST88_T)[3:ncol(rela_ST88_T)] = newcol
rela_ST88_T = rela_ST88_T[,!colnames(rela_ST88_T) %in% NA]
col_rmNA = colnames(rela_ST88_T)[3:ncol(rela_ST88_T)]
#step4：取均值ֵ
avg_ST88_T = data.frame(Index = col_rmNA, avg = apply(rela_ST88_T[,3:ncol(rela_ST88_T)], 2, mean))



########################################################################################################################
library(fs)
#Set the path for storing output files.
dir.create('share_metabolic/')
setwd('share_metabolic/')

#output1 = paste0(sample_name,'.rela.csv')
#write.csv(rela_ST88, output1, row.names=F)
#output2 = paste0(sample_name,'.rela_T.csv')
#write.csv(rela_ST88_T, output2, row.names=F)
output3 = paste0(sample_name,'.avg_T.csv')
write.csv(avg_ST88_T, output3, row.names=F)













