
#Visualization of density distribution of m/zs.


#You can get a list of 'plot.kde.index_{index}.txt' from 'TANG.merge.mass.pl'.
#Create a new file named 'list.txt' that each line is the name of 'plot.kde.index_{index}.txt', such as 'plot.kde.index_1001.txt'.

file_list <- read.table(file="/data/bingling/rPackageTutorial/Continued/example/Dataset/list.txt", header=F)$V1

for (i in 1:length(file_list)){
    output_name = sub('.txt', '.pdf', file_list[i])
    tab <- as.data.frame(read.table(file=paste0(file_list[i]), header=F, sep="\t"))

    KDE <- subset(tab, V3=="KDE")
    Cluster <- subset(tab, V3=="Cluster")
    Mass <- subset(tab, V3=="Mass")

    pdf(file=output_name, width=10, height=6)
    p1 <- plot(KDE$V1, KDE$V2, type='h')
    p2 <- points(Cluster$V1, Cluster$V2, col="red", cex=1)
    p3 <- points(Mass$V1, Mass$V2, col="blue", cex=1)

    print(p1)
    print(p2)
    print(p3)

    dev.off()
}
