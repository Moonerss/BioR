#' PCA analysis of data
#'
#' @importFrom ggpubr stat_conf_ellipse
#' @param data A matrix representing the genomic data such as gene expression data, miRNA expression data.\cr
#' For the matrix, the rows represent the genomic features, and the columns represent the samples.
#' @param group A data frame contain two columns. The first column is sample name matched with colnames of data,
#' The second column is the cluster label of samples.
#' @param pic_title The title of plot.
#' @param ellipse Whether add the confidence ellipse.
#'
#' @return
#' Return a `ggplot` object contained a PCA plot.
#'
#' @export
#'
plot_pca <- function(data, group, pic_title = "All-PCA", ellipse = FALSE) {
  # check samples
  colnames(group) <- c('ID', 'Type')
  if (ncol(data) != length(unique(group$ID))) {
    stop('The sample in `group` not all matched with `data`')
  }

  ##
  ID <- unique(as.vector(group$ID))
  subdata <- subset(data, select = ID)
  #subdata=log2(subdata)  #定量矩阵有0值，不能直接log转换
  subdata <- t(na.omit(subdata))
  data.pca <- prcomp(subdata, scale. = TRUE)
  # summary(data.pca)
  pca <- as.data.frame(data.pca$x)
  pca$ID <- rownames(pca)
  pca <- merge(pca, group, by = "ID")
  percent1 <- sprintf("%.1f%%", data.pca$sdev[1]^2/sum(data.pca$sdev^2)*100)
  percent2 <- sprintf("%.1f%%", data.pca$sdev[2]^2/sum(data.pca$sdev^2)*100)
  xlabel <- paste("PC1 (",percent1,")")
  ylabel <- paste("PC2 (",percent2,")")

  p <- ggplot(pca) +
	aes(x = PC1, y = PC2, fill = Type) +
	geom_point(shape = "circle filled", size = 2.5, colour = 'black') +
	# geom_text(aes(label = ID), hjust = 0, vjust = 0) +
	scale_fill_brewer(palette = "Set1", direction = 1) +
	labs(x = xlabel, y = ylabel, title = pic_title) +
	guides(fill = guide_legend(title = "Sample type")) +
	theme_bw() +
	theme(plot.title = element_text(hjust=0.5)) 
  
  if (ellipse) {
    p <- p +
      stat_ellipse(aes(colour = .data$Type), level = 0.95, geom = "polygon",
                   alpha = 0.3, lwd = 1, show.legend = F) +
      scale_color_brewer(palette = "Set1", direction = 1)
    # add confidence ellipses use stat_conf_ellipse from ggpubr
      # stat_conf_ellipse(alpha = 0.3, geom = 'polygon')
  }
  return(p)
}
