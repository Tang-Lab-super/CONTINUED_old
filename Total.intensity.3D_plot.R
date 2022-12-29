

library("plotly")
library("viridis")
library("RColorBrewer")


sample_list <- read.table(file="sample.list.txt", header=F)$V1

for (sample_name in sample_list){
  
    setwd(paste0(csv_output_dir_path, sample_name))
    getwd()
    print(sample_name)
########################################################################################################################
# read files
    input_unlockmass <- paste0(sample_name, "_3k_signal.unlock_mass.txt")
    signal_total_unlock <- read.table(input_unlockmass, header=F)
    sig_mat_df <- data.frame(X=as.numeric(signal_total_unlock[2:length(signal_total_unlock$V2), 2]), Y=as.numeric(signal_total_unlock[2:length(signal_total_unlock$V3), 3]))
    sig_mat_df$sum_intensity <- rowSums(signal_total_unlock[2:nrow(signal_total_unlock), 4:ncol(signal_total_unlock)])

########################################################################################################################
    sig_mat_df %>% tidyr::spread(key=Y, value=sum_intensity) -> sig_mat
    sig_mat <- sig_mat[3:nrow(sig_mat)-1,3:ncol(sig_mat)-1]

    sig_mat <- as.matrix(sig_mat)
    colnames(sig_mat) <- 1:dim(sig_mat)[2]
    rownames(sig_mat) <- 1:dim(sig_mat)[1]

###########################################
    dir <- paste0(output_dir, '/', sample_name, '/3D_intensity')
    dir.create(dir)

    fig <- plot_ly(x=1:dim(sig_mat)[1], y=1:dim(sig_mat)[2], z=sig_mat, colors=rev(RColorBrewer::brewer.pal(11, "Spectral"))) %>% add_surface()
    htmlwidgets::saveWidget(fig, file = paste0(dir, '/', sample_name, ".v1.html"), selfcontained =F)

    fig <- plot_ly(x=1:dim(sig_mat)[1], y=1:dim(sig_mat)[2], z=sig_mat) %>% add_surface()
    htmlwidgets::saveWidget(fig, file = paste0(dir, '/', sample_name, ".v2.html"), selfcontained =F)

########################################################################################################################

    sig_kde2d <-list(x=1:dim(sig_mat)[1], y=1:dim(sig_mat)[2], z=sig_mat)

    output_pdf <- paste0(dir, '/', sample_name, "contour.plot.v1.pdf")
    pdf(output_pdf, width=10 * dim(sig_mat)[1] / dim(sig_mat)[2] *1.055, height=10 * dim(sig_mat)[1] / dim(sig_mat)[2])
    cols <- brewer.pal(11, "Spectral")
    cols <- rev(colorRampPalette(cols)(28))
    p1 <- filled.contour(sig_kde2d, col=cols, plot.axes={contour(sig_kde2d, add=T, lwd=1)})
    print(p1)
    dev.off()

    output_pdf <- paste0(dir, '/', sample_name, "contour.plot.v2.pdf")
    pdf(output_pdf, width=10 * dim(sig_mat)[1] / dim(sig_mat)[2] *1.055, height=10 * dim(sig_mat)[1] / dim(sig_mat)[2])
    cols <- viridis_pal(alpha = 1, begin = 0, end = 1, direction = 1, option = "D")(11)
    cols <- colorRampPalette(cols)(28)
    p2 <- filled.contour(sig_kde2d, col=cols, plot.axes={contour(sig_kde2d, add=T, lwd=1)})
    print(p2)
    dev.off()

}












