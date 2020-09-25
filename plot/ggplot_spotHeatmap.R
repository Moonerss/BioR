#' @@ ggplot_spotHeatmap
#' 绘制一个点组成的矩阵，点的大小，颜色，形状都可以映射相应参数
#' @param  df 数据框
#' @param  x_idx 数据框列名，映射x轴
#' @param  y_idx 数据框列名，映射y轴
#' @param  size_idx 数据框列名，映射点的大小
#' @param  fill_idx  数据框列名，映射点的填充颜色，
#' @param  defined_fill 自定义填充颜色，默认fill_idx列是连续型数值
#' @param  is.fill_discrete fill_idx列是否是离散型数值
#' @param  color_idx 数据框列名，映射点的边框颜色
#' @param  defined_color 自定义点的边框颜色
#' @param  Angle 坐标轴是否倾斜
#' @param  output_file 绘图结果存储文件
#' @param  width, height 指定输出PDF的宽和高
#'
#' @returnType 
#' @return 
#' 
#' @author 
#' 
ggplot_spotHeatmap <- function(df, x_idx = "", y_idx = "", color_idx = NULL, fill_idx = "", 
                               defined_fill = list(low = "steelblue", mid = "white", high = "red"),
                               is.fill_discrete = FALSE,
                               defined_color = c("#00AFBB", "#FC4E07"),  
                               size_idx = NULL, alpha_idx = NULL,
                               axis_angle = TRUE,
                               title_text = "", x_lab = "", y_lab = "", 
                               output_file, width, height){
  library(ggplot2)
  library(ggpubr)
  theme_set(theme_pubclean())
  
  ###（1）定义图像主体，映射xy轴，点的填充颜色
  if(is.fill_discrete){
    p_fill <- scale_fill_discrete(defined_fill)
  }else{
    p_fill <- scale_fill_gradientn(colours = defined_fill)
  }
  
  p_color <- scale_colour_manual(values = defined_color)
  #如果定义了边框颜色
  if(!is.null(color_idx)){
    p_main <- ggplot(df, aes(x = df[, x_idx], y = df[,y_idx], color = df[, color_idx], fill = df[, fill_idx]))
    p_assembly <- p_main + p_fill + p_color
  }else{
    p_main <- ggplot(df, aes(x = df[, x_idx], y = df[,y_idx], color = NULL, fill = df[, fill_idx]))
    p_assembly <- p_main + p_fill
  }
  
  ###（2）点的大小和透明度映射一列特征，也可直接指定（例如所有点大小为6）
  if(is.null(size_idx) & is.null(alpha_idx)){
    p_point <- geom_point(shape = 21, size = 6)
  }else if(is.null(alpha_idx)){p_point <- geom_point(shape = 21, aes(size = df[,size_idx]))
  }else{
    p_point <- geom_point(shape = 21, aes(size = df[,size_idx], alpha = df[,alpha_idx])) 
  }
  p_assembly <- p_assembly + p_point
  
  ###（3）坐标轴和坐标名称
  p_axis <- theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1))
  p_labs <- labs(x = x_lab, y = y_lab, title = title_text)
  #如果倾斜x轴标签
  if(axis_angle){
    p_assembly <- p_assembly + p_axis + p_labs
  }else{
    p_assembly <- p_assembly + p_labs
  }
  
  ###（4）legend label
  p_assembly <- p_assembly + guides(size = guide_legend(size_idx), fill = guide_colourbar(fill_idx))
  
  ggsave(filename = output_file, p_assembly, width = width, height = height)
  return(p_assembly)
 
}