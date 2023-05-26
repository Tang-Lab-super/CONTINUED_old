
library(WGCNA)


# Make a WGCNA input matrix for the tumor region.
lT=c('ST103_20210718','ST109_20210330','ST121_20210806','ST124_20211223',
'ST129_20201214','ST129_20210428','ST133_20210429','ST32_20210807','ST35_20210401',
'ST49_20211220','ST69_20211222','ST84_20211223','ST87_20210331','ST88_20210331')
sample_name = 'ST06_20210716'
file_dir = paste0('share_metabolic/',sample_name,'.avg_T.csv')
inten_tumor = read.csv(file_dir)
colnames(inten_tumor)[2] = paste0(sample_name,'_avg')
intensity_tumor_merge = inten_tumor  #构造初始
for (sample_name in lT){
  print(sample_name)
  file_dir = paste0('share_metabolic/',sample_name,'.avg_T.csv')
  inten_tumor = read.csv(file_dir)
  colnames(inten_tumor)[2] = paste0(sample_name,'_avg')
  intensity_tumor_merge = merge(intensity_tumor_merge,inten_tumor,by.x='Index',by.y='Index',all=TRUE)
  print(dim(intensity_tumor_merge))
}
rownames(intensity_tumor_merge) = intensity_tumor_merge$Index
intensity_tumor_merge = intensity_tumor_merge[,-1]
intensity_tumor_merge[is.na(intensity_tumor_merge)]=0
intensity_tumor_merge = t(intensity_tumor_merge)
mat_expr_selected = intensity_tumor_merge

#######################################################################
# Make a WGCNA input matrix for the tumor border region.
lTB=c('ST103_20210718','ST109_20210330','ST121_20210806','ST124_20211223',
'ST129_20201214','ST129_20210428','ST133_20210429','ST32_20210807',
'ST49_20211220','ST69_20211222','ST87_20210331')
sample_name = 'ST88_20210331'
file_dir = paste0('share_metabolic/',sample_name,'.avg_TB.csv')
inten_tumor = read.csv(file_dir)
colnames(inten_tumor)[2] = paste0(sample_name,'_avg')
intensity_tumor_merge = inten_tumor  #构造初始
for (sample_name in lTB){
  print(sample_name)
  file_dir = paste0('share_metabolic/',sample_name,'.avg_TB.csv')
  inten_tumor = read.csv(file_dir)
  colnames(inten_tumor)[2] = paste0(sample_name,'_avg')
  intensity_tumor_merge = merge(intensity_tumor_merge,inten_tumor,by.x='Index',by.y='Index',all=TRUE)
  print(dim(intensity_tumor_merge))
}
rownames(intensity_tumor_merge) = intensity_tumor_merge$Index
intensity_tumor_merge = intensity_tumor_merge[,-1]
intensity_tumor_merge[is.na(intensity_tumor_merge)]=0
intensity_tumor_merge = t(intensity_tumor_merge)
mat_expr_selected = intensity_tumor_merge



#######################################################################
# Make a WGCNA input matrix for the mucosa region.
lmu=c('ST103_20210718','ST109_20210330','ST121_20210806','ST124_20211223',
'ST129_20201214','ST129_20210428','ST133_20210429','ST32_20210807','ST35_20210401',
'ST49_20211220','ST69_20211222','ST84_20211223','ST87_20210331','ST88_20210331')
sample_name = 'ST06_20210716'
file_dir = paste0('share_metabolic/',sample_name,'.avg_mu.csv')
inten_tumor = read.csv(file_dir)
colnames(inten_tumor)[2] = paste0(sample_name,'_avg')
intensity_tumor_merge = inten_tumor  #构造初始
for (sample_name in lmu){
  print(sample_name)
  file_dir = paste0('share_metabolic/',sample_name,'.avg_mu.csv')
  inten_tumor = read.csv(file_dir)
  colnames(inten_tumor)[2] = paste0(sample_name,'_avg')
  intensity_tumor_merge = merge(intensity_tumor_merge,inten_tumor,by.x='Index',by.y='Index',all=TRUE)
  print(dim(intensity_tumor_merge))
}
rownames(intensity_tumor_merge) = intensity_tumor_merge$Index
intensity_tumor_merge = intensity_tumor_merge[,-1]
intensity_tumor_merge[is.na(intensity_tumor_merge)]=0
intensity_tumor_merge = t(intensity_tumor_merge)
mat_expr_selected = intensity_tumor_merge

########################### WGCNA ###########################

## power
powers <- 1:50
sft <- pickSoftThreshold(mat_expr_selected, powerVector = powers, verbose = 5)


pdf("WGCNA.fitIndices.pdf",width = 12)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], type = 'n',
    xlab = 'Soft Threshold (power)', ylab = 'Scale Free Topology Model Fit, signed R^2',
    main = paste('Scale independence'))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels = powers, col = 'red');
abline(h = 0.90, col = 'red')
dev.off()

pdf("WGCNA.connectivity.pdf")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab = 'Soft Threshold (power)', ylab = 'Mean Connectivity', type = 'n',
    main = paste('Mean connectivity'))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = 'red')
dev.off()


head(sft$fitIndices, n=10)




#上一步估计的最佳 power 值
powers <- sft$powerEstimate
#若无最佳powers值，暂时取经验值9
powers = 9
#检验选定的β值下记忆网络是否逼近 scale free
k <- softConnectivity(datE=mat_expr_selected,power=powers)
sizeGrWindow(10, 5)
pdf("WGCNA.scaleFreePlot.pdf",width=12)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()

#获得 TOM 矩阵
adjacency <- adjacency(mat_expr_selected, power = powers)
tom_sim <- TOMsimilarity(adjacency)
rownames(tom_sim) <- rownames(adjacency)
colnames(tom_sim) <- colnames(adjacency)
tom_sim[1:6,1:6]



#TOM 相异度 = 1 – TOM 相似度
tom_dis  <- 1 - tom_sim

#层次聚类树，使用中值的非权重成对组法的平均聚合聚类
geneTree <- hclust(as.dist(tom_dis), method = 'average')
pdf("WGCNA.geneTree.pdf")
plot(geneTree, xlab = '', sub = '', main = 'Gene clustering on TOM-based dissimilarity',
    labels = FALSE, hang = 0.04)
dev.off()


#使用动态剪切树挖掘模块
minModuleSize <- 30  #模块基因数目
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = tom_dis,
    deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)

table(dynamicMods)


#模块颜色指代
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

pdf("WGCNA.geneTree.dynamicColors.pdf")
plotDendroAndColors(geneTree, dynamicColors, 'Dynamic Tree Cut',
    dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05,
    main = 'Gene dendrogram and module colors')
dev.off()




#基因表达聚类树和共表达拓扑热图
plot_sim <- -(1-tom_sim)
#plot_sim <- log(tom_sim)
diag(plot_sim) <- NA

pdf("WGCNA.TOMplot.pdf")
TOMplot(plot_sim, geneTree, dynamicColors,
    main = 'Network heatmap plot, selected genes')
dev.off()

net = blockwiseModules(mat_expr_selected, power = powers, maxBlockSize = 5000,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.1,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = "pearson",
                       loadTOMs=TRUE,saveTOMFileBase = "data.tom",
                       verbose = 3)
library(stringr)
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
pdf("WGCNA.plotEigengeneNetworks.pdf")
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)
dev.off()

#合并相似module
#通过dynamicTreeCut识别到modules之后，还会结合每个modules的基因表达量数据，来识别相关性很高的modules, 从而进行合并，其原理是对modules进行聚类
MEDissThres=0.25
merge = mergeCloseModules(mat_expr_selected, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
png("WGCNA.geneTree.dynamicColors.merge.png")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

png("WGCNA.TOMplot.merge.png")
TOMplot(plot_sim, geneTree, mergedColors,
    main = 'Network heatmap plot, selected genes')
dev.off()
table(mergedColors)



############################################################


#统计每个颜色模块的代谢物数量
cl=c()
module_gene=c()
module_colors=unique(mergedColors)
for (color in module_colors){
    module=probes[which(mergedColors==color)]
    group1=rep(color,times=length(module))
    cl=c(cl,group1)
    module_gene=c(module_gene,module)
}
my_module=data.frame("colour"=cl,"gene"=module_gene)
#write.csv(my_module,"Modules_gene.csv")
write.csv(table(mergedColors),'module_colors.csv')


#取出感兴趣模块的代谢物，输入metaboanalyst做代谢通路分析
colon = read.csv(file="/data/bingling/rPackageTutorial/Continued/example/Dataset/colon_cancer_desi.clustered_mass.table.with.anno.csv", header=T)
module_index = my_module[my_module$colour=='lightgreen',]$gene
colon_T3 = colon[colon$Index %in% module_index,]
length(colon_T3[colon_T3$anno_lipid!='None',3])   #[1] 30
length(colon_T3[colon_T3$anno_small_mol!='None',4])   #[1] 21
out_lip3 = get_lipid(colon_T3)
table(gsub('[(].*','',out_lip3[[4]]))
out_mol3 = get_mol(colon_T3)
#https://www.metaboanalyst.ca/MetaboAnalyst/upload/PathUploadView.xhtml




########### WGCNA富集气泡图 ##############
#在metaboanalyst网站上做完通路富集之后，如果想把不同模块的结果拼在一起，则参考运行下面的代码

data_tumor = read.csv('WGCNA.mu.csv')
pdf('WGCNA.mu.dot.pdf',height = 5)
t = c("Arginine biosynthesis", "Citrate cycle (TCA cycle)","Alanine, aspartate and glutamate metabolism","Linoleic acid metabolism",
"D-Glutamine and D-glutamate metabolism","Nitrogen metabolism","Taurine and hypotaurine metabolism","Starch and sucrose metabolism",
"Biosynthesis of unsaturated fatty acids","Butanoate metabolism","Pantothenate and CoA biosynthesis","beta-Alanine metabolism")

t = c("Taurine and hypotaurine metabolism","Sphingolipid metabolism","Glycerophospholipid metabolism","Biosynthesis of unsaturated fatty acids",
"Fatty Acid Biosynthesis","Malate-Aspartate Shuttle","Phosphatidylethanolamine Biosynthesis","Butyrate Metabolism",
"Alpha Linolenic Acid and Linoleic Acid Metabolism")
data_tumor$pathway = factor(data_tumor$pathway, levels=rev(t))
ggplot(data_tumor,aes(x = as.factor(cluster),y = as.factor(pathway)))+
    geom_point(aes(color = pvalue,
                   size = enrichment))+
	scale_size(breaks = seq(20,60,20),limits=c(0,60))+
    scale_color_gradient(low = "red", high = "LemonChiffon",breaks = seq(0,0.15,0.05),limits=c(0,0.15),na.value='LemonChiffon')+
	#xlab("Fold Enrichment")+
    theme_bw()+
	theme(panel.grid=element_blank())+
	xlab(NULL)+
	ylab(NULL)
	#guides(color=guide_legend(order = 1),size=guide_legend(order=0))
    #guides(
        #reverse color order (higher value on top)
        #color = guide_colorbar(reverse = TRUE))
        #reverse size order (higher diameter on top)
        #size = guide_legend(reverse = TRUE))
dev.off()

pdf('WGCNA_mu.dot.pdf',height = 5)
ggplot(data_tumor,aes(x = as.factor(cluster),y = pathway))+
    geom_point(aes(color = pvalue,
                   size = enrichment))+
    scale_color_gradient(low = "red", high = "LemonChiffon")+
	#xlab("Fold Enrichment")+
    theme_bw()+
	theme(panel.grid=element_blank())+
    #guides(
        #reverse color order (higher value on top)
        #color = guide_colorbar(reverse = TRUE))
        #reverse size order (higher diameter on top)
        #size = guide_legend(reverse = TRUE))

dev.off()










