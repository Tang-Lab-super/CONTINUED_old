
#Draw a three-dimensional map of the metabolites of the omega6 metabolic pathway: LA, DGLA, AA, ADA

library(plotly)
library(tidyverse)
library(gmodels)
library(plot3D)

infor <- as.data.frame(read.table(file="file.infor.txt", header=T))
dir_path <- '/data/tang/desi/processing/supplementary/'   #obj_path的上级目录
dir_mass_merged <- '/data/tang/desi/processing/merge.lockmass/merged.mass.for.sample.cytof/'   #'{sample_name}_mass.txt'的上级目录

####################### ST129 ########################
# get matrix from sample 'ST129_20210428'
l=17
sample_name <- infor[l, 1]
obj_path <- infor[l, 2]
assign(paste0('obj_', sample_name), Create_New_Obj(sample_name, dir_path, obj_path, dir_mass_merged))

input_unlockmass <- paste0(sample_name, "_signal.unlock_mass.txt")
signal_total_unlock <- read.table(input_unlockmass, header=F)
signal_total_unlock[2:length(signal_total_unlock$V2),]$V2 <- as.character(as.numeric(signal_total_unlock[2:length(signal_total_unlock$V2),]$V2))
signal_total_unlock[2:length(signal_total_unlock$V3),]$V3 <- as.character(as.numeric(signal_total_unlock[2:length(signal_total_unlock$V3),]$V3))
sig_mat_df <- data.frame(X=as.numeric(signal_total_unlock[2:length(signal_total_unlock$V2), 2]), Y=as.numeric(signal_total_unlock[2:length(signal_total_unlock$V3), 3]))

dir_name <- paste0('./', sample_name)
mass_tag <- 'Mass_clustered'
coor = read.csv(file='ST129_20210428_signal.selected.coor.txt', sep='\t')

#准备绘图数据框

#LA:279
grep('279', as.numeric(signal_total_unlock[1,4:ncol(signal_total_unlock)]),value = T)  #[1] "279.218"  "279.0252"
which(as.numeric(signal_total_unlock[1,])=="279.218")  #[1] 8
sig_mat_df$sum_intensity = signal_total_unlock$V8[2:nrow(signal_total_unlock)]
#筛选组织坐标
sig_mat_df_279new = get_tissue_coor(sig_mat_df, coor)   #这一步会运行比较久
#与其他代谢物3D图重合时，若信号丰度差距太大，则多运行下面一行代码，将丰度值压缩到log
#sig_mat_df_279new[sig_mat_df_279new$sum_intensity!=0,]$sum_intensity = log2(sig_mat_df_279new[sig_mat_df_279new$sum_intensity!=0,]$sum_intensity)
sig_mat_df_279new %>% tidyr::spread(key=Y, value=sum_intensity) -> sig_mat
sig_mat <- sig_mat[3:nrow(sig_mat)-1,3:ncol(sig_mat)-1]
sig_mat <- as.matrix(sig_mat)
colnames(sig_mat) <- 1:dim(sig_mat)[2]
rownames(sig_mat) <- 1:dim(sig_mat)[1]
#将生成的文件夹同时下载，网页显示3D动图
dir = '3D/'
sig_mat = t(sig_mat)
fig <- plot_ly(x=1:dim(sig_mat)[1], y=1:dim(sig_mat)[2], z=sig_mat, colors=c('white',brewer.pal(9, "OrRd"))) %>% add_surface()
htmlwidgets::saveWidget(fig, file = paste0(dir, '/', sample_name, ".v1.279.html"), selfcontained =F)


#DGLA:305
grep('305', as.numeric(signal_total_unlock[1,4:ncol(signal_total_unlock)]),value = T)  #[1] "305.2311" "824.5305"
which(as.numeric(signal_total_unlock[1,])=="305.2311")  #[1] 8
sig_mat_df$sum_intensity = signal_total_unlock$V36[2:nrow(signal_total_unlock)]
#筛选组织坐标
sig_mat_df_305new = get_tissue_coor(sig_mat_df, coor)
sig_mat_df_305new %>% tidyr::spread(key=Y, value=sum_intensity) -> sig_mat1
sig_mat1 <- sig_mat1[3:nrow(sig_mat1)-1,3:ncol(sig_mat1)-1]
sig_mat1 <- as.matrix(sig_mat1)
colnames(sig_mat1) <- 1:dim(sig_mat1)[2]
rownames(sig_mat1) <- 1:dim(sig_mat1)[1]

#AA:303
grep('303', as.numeric(signal_total_unlock[1,4:ncol(signal_total_unlock)]),value = T)  #[1] "303.2162" "962.4303"
which(as.numeric(signal_total_unlock[1,])=="303.2162")  #[1] 12
sig_mat_df$sum_intensity = signal_total_unlock$V12[2:nrow(signal_total_unlock)]
sig_mat_df_303new = get_tissue_coor(sig_mat_df, coor)
sig_mat_df_303new %>% tidyr::spread(key=Y, value=sum_intensity) -> sig_mat2
sig_mat2 <- sig_mat2[3:nrow(sig_mat2)-1,3:ncol(sig_mat2)-1]
sig_mat2 <- as.matrix(sig_mat2)
colnames(sig_mat2) <- 1:dim(sig_mat2)[2]
rownames(sig_mat2) <- 1:dim(sig_mat2)[1]

#ADA:331
grep('331', as.numeric(signal_total_unlock[1,4:ncol(signal_total_unlock)]),value = T)  #[1] "331.2456" "331.1445"
which(as.numeric(signal_total_unlock[1,])=="331.2456")  #[1] 42
sig_mat_df$sum_intensity = signal_total_unlock$V42[2:nrow(signal_total_unlock)]
sig_mat_df_331new = get_tissue_coor(sig_mat_df, coor0)
sig_mat_df_331new %>% tidyr::spread(key=Y, value=sum_intensity) -> sig_mat3
sig_mat3 <- sig_mat3[3:nrow(sig_mat3)-1,3:ncol(sig_mat3)-1]
sig_mat3 <- as.matrix(sig_mat3)
colnames(sig_mat3) <- 1:dim(sig_mat3)[2]
rownames(sig_mat3) <- 1:dim(sig_mat3)[1]


sig_mat[sig_mat==0]=NA
sig_mat1[sig_mat1==0]=NA
sig_mat2[sig_mat2==0]=NA
#只保留ADA丰度最高的坐标点
sig_mat3[sig_mat3<2300]=0
sig_mat3[sig_mat3==0]=NA
library(plot3D)
library(plotly)

#单独绘制
pdf('3D/ST129.LA.pdf',width=10)
persp3D(z = sig_mat,  clim = c(0,9000),col = ramp.col(col = c('white','Dodgerblue')),box=FALSE, theta = 0, phi = 90)
dev.off()
pdf('3D/ST129.DGLA.pdf',width=10)
persp3D(z = sig_mat1, col = ramp.col(col = c('white','Tomato')),box=FALSE, theta = 0, phi = 90)
dev.off()
pdf('3D/ST129.AA.pdf',width=10)
persp3D(z = sig_mat2, col = ramp.col(col = c('white','DarkOrange')),box=FALSE,theta = 0, phi = 90)
dev.off()
pdf('3D/ST129.ADA.pdf',width=10)
persp3D(z = sig_mat3, col = ramp.col(col = c('white','Magenta1')),box=FALSE,theta = 0, phi = 90)
dev.off()

#绘制LA和ADA重合的3D图
pdf('3D/ST129.test2.pdf',width=10)
persp3D(z = sig_mat, clim = c(0,9000),col = ramp.col(col = c('white','Dodgerblue')),box=FALSE, plot =FALSE, theta = 0, phi = 90)
persp3D(z = sig_mat3, clim = c(0,9000),col = ramp.col(col = c('white','Magenta1')),box=FALSE, theta = 0, phi = 90,add = TRUE)
dev.off()


##############################ST103
l=29
sample_name <- infor[l, 1]
obj_path <- infor[l, 2]
assign(paste0('obj_', sample_name), Create_New_Obj(sample_name, dir_path, obj_path, dir_maass_merged))

input_unlockmass <- paste0(sample_name, "_signal.unlock_mass.txt")
signal_total_unlock <- read.table(input_unlockmass, header=F)
signal_total_unlock[2:length(signal_total_unlock$V2),]$V2 <- as.character(as.numeric(signal_total_unlock[2:length(signal_total_unlock$V2),]$V2))
signal_total_unlock[2:length(signal_total_unlock$V3),]$V3 <- as.character(as.numeric(signal_total_unlock[2:length(signal_total_unlock$V3),]$V3))
sig_mat_df <- data.frame(X=as.numeric(signal_total_unlock[2:length(signal_total_unlock$V2), 2]), Y=as.numeric(signal_total_unlock[2:length(signal_total_unlock$V3), 3]))
coor = read.csv(file='ST103_20210718_signal.selected.coor.txt', sep='\t')

#LA:279
grep('279', as.numeric(signal_total_unlock[1,4:ncol(signal_total_unlock)]),value = T)  #[1] "279.2886" "279.0946" "869.2796"
which(as.numeric(signal_total_unlock[1,])=="279.2886")  #[1] 9
sig_mat_df$sum_intensity = signal_total_unlock$V9[2:nrow(signal_total_unlock)]
#筛选组织坐标
sig_mat_df_279_ST103 = get_tissue_coor(sig_mat_df, coor)
#sig_mat_df_279_ST103$sum_intensity = sig_mat_df_279_ST103$sum_intensity/5
sig_mat_df_279_ST103 %>% tidyr::spread(key=Y, value=sum_intensity) -> sig_mat
sig_mat <- sig_mat[3:nrow(sig_mat)-1,3:ncol(sig_mat)-1]
sig_mat <- as.matrix(sig_mat)
colnames(sig_mat) <- 1:dim(sig_mat)[2]
rownames(sig_mat) <- 1:dim(sig_mat)[1]

#DGLA:305
grep('305', as.numeric(signal_total_unlock[1,4:ncol(signal_total_unlock)]),value = T)  #[1] "305.3059" "504.3057"
which(as.numeric(signal_total_unlock[1,])=="305.3059")  #[1] 315
sig_mat_df$sum_intensity = signal_total_unlock$V315[2:nrow(signal_total_unlock)]
#筛选组织坐标
sig_mat_df_305new = get_tissue_coor(sig_mat_df, coor)
sig_mat_df_305new %>% tidyr::spread(key=Y, value=sum_intensity) -> sig_mat1
sig_mat1 <- sig_mat1[3:nrow(sig_mat1)-1,3:ncol(sig_mat1)-1]
sig_mat1 <- as.matrix(sig_mat1)
colnames(sig_mat1) <- 1:dim(sig_mat1)[2]
rownames(sig_mat1) <- 1:dim(sig_mat1)[1]


#AA:303
grep('303', as.numeric(signal_total_unlock[1,4:ncol(signal_total_unlock)]),value = T)  #[1] "303.291"  "303.1083"
which(as.numeric(signal_total_unlock[1,])=="303.291")  #[1] 47
sig_mat_df$sum_intensity = signal_total_unlock$V47[2:nrow(signal_total_unlock)]
sig_mat_df_303new = get_tissue_coor(sig_mat_df, coor)
sig_mat_df_303new %>% tidyr::spread(key=Y, value=sum_intensity) -> sig_mat2
sig_mat2 <- sig_mat2[3:nrow(sig_mat2)-1,3:ncol(sig_mat2)-1]
sig_mat2 <- as.matrix(sig_mat2)
colnames(sig_mat2) <- 1:dim(sig_mat2)[2]
rownames(sig_mat2) <- 1:dim(sig_mat2)[1]


#ADA:331
grep('331', as.numeric(signal_total_unlock[1,4:ncol(signal_total_unlock)]),value = T)  #[1] "1040.9331" "331.3246"
which(as.numeric(signal_total_unlock[1,])=="331.3246")  #[1] 399
sig_mat_df$sum_intensity = signal_total_unlock$V399[2:nrow(signal_total_unlock)]
sig_mat_df_331new = get_tissue_coor(sig_mat_df, coor)
sig_mat_df_331new %>% tidyr::spread(key=Y, value=sum_intensity) -> sig_mat3
sig_mat3 <- sig_mat3[3:nrow(sig_mat3)-1,3:ncol(sig_mat3)-1]
sig_mat3 <- as.matrix(sig_mat3)
colnames(sig_mat3) <- 1:dim(sig_mat3)[2]
rownames(sig_mat3) <- 1:dim(sig_mat3)[1]
#绘制3D动态图，确定阈值600
dir = '3D/'
sig_mat3 = t(sig_mat3)
fig <- plot_ly(x=1:dim(sig_mat3)[1], y=1:dim(sig_mat3)[2], z=sig_mat3, colors=c('white',brewer.pal(9, "OrRd"))) %>% add_surface()
htmlwidgets::saveWidget(fig, file = paste0(dir, '/', sample_name, ".v1.331.html"), selfcontained =F)
sig_mat3[sig_mat3<600]=0


sig_mat[sig_mat==0]=NA
sig_mat1[sig_mat1==0]=NA
sig_mat2[sig_mat2==0]=NA
sig_mat3[sig_mat3<600]=0
sig_mat3[sig_mat3==0]=NA
library(plot3D)
library(plotly)
#分别绘制
pdf('3D/ST103.LA.pdf',width=10)
persp3D(z = sig_mat, clim = c(0,10000),col = ramp.col(col = c('white','Dodgerblue')),box=FALSE, theta = 0, phi = 90)
dev.off()
pdf('3D/ST103.DGLA.pdf',width=10)
persp3D(z = sig_mat1, col = ramp.col(col = c('white','Tomato')),box=FALSE, theta = 0, phi = 90)
dev.off()
pdf('3D/ST103.AA.pdf',width=10)
persp3D(z = sig_mat2, col = ramp.col(col = c('white','DarkOrange')),box=FALSE,theta = 0, phi = 90)
dev.off()
pdf('3D/ST103.ADA.pdf',width=10)
persp3D(z = sig_mat3, col = ramp.col(col = c('white','Magenta1')),box=FALSE,theta = 0, phi = 90)
dev.off()
#绘制LA和ADA重合的3D图
pdf('3D/ST103.pdf',width=10)
persp3D(z = sig_mat, clim = c(0,2000),col = ramp.col(col = c('white','Dodgerblue')),box=FALSE, plot =FALSE, theta = 0, phi = 90)
persp3D(z = sig_mat3, clim = c(0,2000),col = ramp.col(col = c('white','Magenta1')),box=FALSE, theta = 0, phi = 90,add = TRUE)
dev.off()


##############################ST32
l=20
sample_name <- infor[l, 1]
obj_path <- infor[l, 2]
assign(paste0('obj_', sample_name), Create_New_Obj(sample_name, dir_path, obj_path, dir_maass_merged))

input_unlockmass <- paste0(sample_name, "_signal.unlock_mass.txt")
signal_total_unlock <- read.table(input_unlockmass, header=F)
signal_total_unlock[2:length(signal_total_unlock$V2),]$V2 <- as.character(as.numeric(signal_total_unlock[2:length(signal_total_unlock$V2),]$V2))
signal_total_unlock[2:length(signal_total_unlock$V3),]$V3 <- as.character(as.numeric(signal_total_unlock[2:length(signal_total_unlock$V3),]$V3))
sig_mat_df <- data.frame(X=as.numeric(signal_total_unlock[2:length(signal_total_unlock$V2), 2]), Y=as.numeric(signal_total_unlock[2:length(signal_total_unlock$V3), 3]))
coor = read.csv(file='ST32_20210807_signal.selected.coor.txt', sep='\t')


#LA
grep('279', as.numeric(signal_total_unlock[1,4:ncol(signal_total_unlock)]),value = T)  #[1] "279.3069" "279.1198" "635.6279" "863.6279"
which(as.numeric(signal_total_unlock[1,])=="279.3069")  #[1] 7
sig_mat_df$sum_intensity = signal_total_unlock$V7[2:nrow(signal_total_unlock)]
sig_mat_df_279_ST103 = get_tissue_coor(sig_mat_df, coor)
#sig_mat_df_279_ST103$sum_intensity = sig_mat_df_279_ST103$sum_intensity/3.5
sig_mat_df_279_ST103 %>% tidyr::spread(key=Y, value=sum_intensity) -> sig_mat
sig_mat <- sig_mat[3:nrow(sig_mat)-1,3:ncol(sig_mat)-1]
sig_mat <- as.matrix(sig_mat)
colnames(sig_mat) <- 1:dim(sig_mat)[2]
rownames(sig_mat) <- 1:dim(sig_mat)[1]

#DGLA
grep('305', as.numeric(signal_total_unlock[1,4:ncol(signal_total_unlock)]),value = T)  #[1] "305.3059" "504.3057"
which(as.numeric(signal_total_unlock[1,])=="305.3252")  #[1] 134
sig_mat_df$sum_intensity = signal_total_unlock$V134[2:nrow(signal_total_unlock)]
sig_mat_df_305new = get_tissue_coor(sig_mat_df, coor)
sig_mat_df_305new %>% tidyr::spread(key=Y, value=sum_intensity) -> sig_mat1
sig_mat1 <- sig_mat1[3:nrow(sig_mat1)-1,3:ncol(sig_mat1)-1]
sig_mat1 <- as.matrix(sig_mat1)
colnames(sig_mat1) <- 1:dim(sig_mat1)[2]
rownames(sig_mat1) <- 1:dim(sig_mat1)[1]

#AA
grep('303', as.numeric(signal_total_unlock[1,4:ncol(signal_total_unlock)]),value = T)  #[1] "303.3103" "295.3036" "303.1486"
which(as.numeric(signal_total_unlock[1,])=="303.3103")  #[1] 13
sig_mat_df$sum_intensity = signal_total_unlock$V13[2:nrow(signal_total_unlock)]
sig_mat_df_303new = get_tissue_coor(sig_mat_df, coor)
sig_mat_df_303new %>% tidyr::spread(key=Y, value=sum_intensity) -> sig_mat2
sig_mat2 <- sig_mat2[3:nrow(sig_mat2)-1,3:ncol(sig_mat2)-1]
sig_mat2 <- as.matrix(sig_mat2)
colnames(sig_mat2) <- 1:dim(sig_mat2)[2]
rownames(sig_mat2) <- 1:dim(sig_mat2)[1]

#ADA
grep('331', as.numeric(signal_total_unlock[1,4:ncol(signal_total_unlock)]),value = T)  #[1] "389.3331" "331.3452" "800.6331" "329.3311"
which(as.numeric(signal_total_unlock[1,])=="331.3452")  #[1] 80
sig_mat_df$sum_intensity = signal_total_unlock$V80[2:nrow(signal_total_unlock)]
sig_mat_df_331new = get_tissue_coor(sig_mat_df, coor)
sig_mat_df_331new %>% tidyr::spread(key=Y, value=sum_intensity) -> sig_mat3
sig_mat3 <- sig_mat3[3:nrow(sig_mat3)-1,3:ncol(sig_mat3)-1]
sig_mat3 <- as.matrix(sig_mat3)
colnames(sig_mat3) <- 1:dim(sig_mat3)[2]
rownames(sig_mat3) <- 1:dim(sig_mat3)[1]
#绘制3D动态图，确定阈值5000
dir = '3D/'
sig_mat3 = t(sig_mat3)
fig <- plot_ly(x=1:dim(sig_mat3)[1], y=1:dim(sig_mat3)[2], z=sig_mat3, colors=c('white',brewer.pal(9, "OrRd"))) %>% add_surface()
htmlwidgets::saveWidget(fig, file = paste0(dir, '/', sample_name, ".v1.331.html"), selfcontained =F)
sig_mat3[sig_mat3<5000]=0


sig_mat[sig_mat==0]=NA
sig_mat1[sig_mat1==0]=NA
sig_mat2[sig_mat2==0]=NA
sig_mat3[sig_mat3<5000]=0
sig_mat3[sig_mat3==0]=NA
library(plot3D)
library(plotly)
#分别绘制
pdf('ST32.LA.pdf',width=10)
persp3D(z = sig_mat, clim = c(0,49000),col = ramp.col(col = c('white','Dodgerblue')),box=FALSE, theta = 0, phi = 90)
dev.off()
pdf('ST32.DGLA.pdf',width=10)
persp3D(z = sig_mat1, col = ramp.col(col = c('white','Tomato')),box=FALSE, theta = 0, phi = 90)
dev.off()
pdf('ST32.AA.pdf',width=10)
persp3D(z = sig_mat2, col = ramp.col(col = c('white','DarkOrange')),box=FALSE,theta = 0, phi = 90)
dev.off()
pdf('ST32.ADA.pdf',width=10)
persp3D(z = sig_mat3, col = ramp.col(col = c('white','Magenta1')),box=FALSE,theta = 0, phi = 90)
dev.off()
#绘制LA和ADA重合的3D图
pdf('3D/ST32.pdf',width=10)
persp3D(z = sig_mat, clim = c(0,14000),col = ramp.col(col = c('white','Dodgerblue')),box=FALSE, plot =FALSE, theta = 0, phi = 90)
persp3D(z = sig_mat3, clim = c(0,14000),col = ramp.col(col = c('white','Magenta1')),box=FALSE, theta = 0, phi = 90,add = TRUE)
dev.off()









