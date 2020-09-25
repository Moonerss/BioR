#' @TODO 利用ggplot输出boxplot_ggpaired
#' @param  inputData_ggplot 输入矩阵
#' @param  outfile 输出文件路径，需要给出扩展名，例如 test.pdf
#' @param  x 使用哪一列作为x坐标轴
#' @param  y 使用哪一列作为y坐标轴
#' @param  group 使用哪一列作为分组依据
#' @param  facet.by 使用哪一列作为分组依据，注意，这个是基于是facet的，设置这个参数后，就会基于facet.by列构建分组面板，然后再设置x即可，保持group=NULL
#' @param  pval_group 基于哪一列作为分组依据，进行统计学检验
#' @param  method 使用什么统计检验，如果设置为NULL则自己判断
#' @param  label 统计学检验在图上的输出模式
# label="p.signif" 输出 *, ***, ns
# label=NULL 输出 Wilcoxon, p=0.005
#' @param  col_palette 颜色模板
#' @param  ylim y轴取值范围
#' @param  fontsize 字体大小
#' @param  width 输出pdf的宽
#' @param  height 输出pdf的高
#' @returnType
#' @return
#'
#' @author ZGX
boxplot_ggpaired <- function(inputData_ggplot=NULL, outfile=NULL, x=NULL, y=NULL, group=NULL, facet.by=NULL, pval_group=NULL, method=c("wilcox.test", "t.test"), label="p.signif", col_palette=c("lancet", "npg", "aaas", "jco", "ucscgb", "uchicago", "simpsons", "rickandmorty"), ylim=NULL, fontsize=16, width=7, height=7){
  library(ggpubr)
  # http://www.sthda.com/english/articles/32-r-graphics-essentials/132-plot-grouped-data-box-plot-bar-plot-and-more/
  ggobj <- ggpaired(inputData_ggplot, x=x, y=y, color=group, palette=col_palette[1], facet.by=facet.by, line.color = "gray", line.size = 0.4) 
  resArr <- NA
  # 是否针对pval_group进行统计学检验
  if(!is.null(pval_group)){
    # ggobj <- ggobj + stat_compare_means(aes(group=pval_group))
    if(!is.null(facet.by)){
      pval_group <- x
    }
    if(is.null(label)){
      drawExprs <- paste0("ggobj <- ggobj + stat_compare_means(aes(paired=TRUE, method='", method[1], "', group=", pval_group, "))")
    }else{
      drawExprs <- paste0("ggobj <- ggobj + stat_compare_means(aes(paired=TRUE, method='", method[1], "', group=", pval_group, "), label='", label, "')")
    }
    # 需要用表达式做。。。
    eval(parse(text=drawExprs))
    if(!is.null(facet.by)){
      eval(parse(text=paste0("resArr <- as.data.frame(compare_means(", y, " ~ ", pval_group, ", data=inputData_ggplot, method='", method[1], "', paired=TRUE, group.by='", facet.by, "'))")))
    }else{
      eval(parse(text=paste0("resArr <- as.data.frame(compare_means(", y, " ~ ", pval_group, ", data=inputData_ggplot, method='", method[1], "', paired=TRUE, group.by='", x, "'))")))
    }
  }
  # 设定主题
  ggobj <- ggobj + theme_bw(base_size=fontsize)
  if(!is.null(ylim))
    ggobj <- ggobj+ylim(ylim)
  # 输出文件
  ggsave(filename=outfile, plot=ggobj, width=width, height=height)
  return(resArr)
}
