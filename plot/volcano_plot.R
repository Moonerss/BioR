#################################################
## TODO:绘制火山图
## Author：Erjie Zhao
## update time： 2020-8-1
#################################################


#' @param diff_expr 差异表达的结果，使用的是limma差异表达的结果
#' @param fdr FDR的取值，默认0.05
#' @param logfc log fold change的取值，默认0.5
#' @param palette 不同散点的颜色
#' @param title 标题名字
#' @param xlab x轴标签
#' @param ylab y轴标签
#' @param mark 是否标记出特殊基因
#' @param mark.logfc 要求logfc大于多少的且fdr小于标准的标记出来
#' @return 返回一个ggplot2对象
#' @import ggplot2, ggrepel
volcano_plot <- function(diff_expr, fdr = 0.05, logfc = 1, palette = c("gray", "red", "darkgreen"),
                         title = NULL, xlab = "log2 fold change", ylab = "-log10 p-value", mark = TRUE, mark.logfc = 1) {
    ## add different labels
    diff_expr <- dplyr::mutate(diff_expr, sig = as.factor(ifelse(diff_expr$adj.P.Val < fdr & abs(diff_expr$logFC) > logfc, 
                                                          ifelse(diff_expr$logFC > logfc,'up','down'),'notsig'))) %>%
                 rownames_to_column(var = "gene")
     ## setting title
     if (is.null(title)) {
          title <- paste0('The number of up gene is ', nrow(diff_expr[diff_expr$sig =='up',]) ,
                          '\nThe number of down gene is ', nrow(diff_expr[diff_expr$sig =='down',]))
     } else {
          title <- title
     }
 
     ## setting color label
     labels <- paste(names(table(diff_expr$sig)), table(diff_expr$sig), sep = " ")
     ## volcano plot
     p <- ggplot(data = diff_expr, aes(x = logFC, y = -log10(adj.P.Val))) +
          ## 将数据分成四象限，对每个象限的数据进行颜色和大小的定义
          geom_point(data = subset(diff_expr, abs(diff_expr$logFC) <= logfc), aes(size = abs(logFC)), color = palette[1], alpha = 0.5) +
          geom_point(data = subset(diff_expr, diff_expr$adj.P.Val >= 0.05 & abs(diff_expr$logFC) > logfc), aes(size = abs(logFC)), color =palette[1], alpha = 0.5) +
          geom_point(data = subset(diff_expr, diff_expr$adj.P.Val < 0.05 & diff_expr$logFC > logfc), aes(size = abs(logFC)), color = palette[2], alpha = 0.5) +
          geom_point(data = subset(diff_expr, diff_expr$adj.P.Val < 0.05 & diff_expr$logFC < -logfc), aes(size = abs(logFC)), color = palette[3], alpha = 0.5) +
          ## 添加分割线（横向+竖向）
          geom_hline(yintercept = -log10(fdr), lty = 4, lwd = 0.6, alpha = 0.8)+
          geom_vline(xintercept = c(-logfc, logfc), lty = 4, lwd = 0.6, alpha = 0.8)+
          #设置背景为空白
          theme_bw()+
          theme(panel.border = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"))+ 
          ## 标题
          labs(title = title, x = xlab, y = ylab, color = "") +
          ## 去掉图例
          theme(legend.position='none')
     ## 是否标记重要基因
     if (mark) {
          p <- p + geom_text_repel(data = dplyr::filter(diff_expr, adj.P.Val < fdr, abs(logFC) > mark.logfc,), aes(label = gene), color = "black", alpha = 0.8)
     }   

    return(p)
}


