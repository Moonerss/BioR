##   TO DO: 转录差异表达
##
##   Author: jeason zhao
#####################################################


#' @name diff_limma
#' @description  使用limma包进行差异表达
#' @param exp.data 表达谱，行为基因，列为样本
#' @param label.list 数据框，一列为样本名称，一列为case和control标签
#' @return 返回limma的差异表达结果

diff_limma <- function(exp.data, label.list) {
  
  ## 载入R包
  suppressPackageStartupMessages(library(limma))
  
  ## 提取分组信息 
  group_list <- label.list[match(label.list[,1], colnames(exp.data)), 2]
  
  ## 构建分组矩阵
  design <- model.matrix(~ 0 + factor(group_list))
  colnames(design) <- levels(factor(group_list))
  rownames(design) <- colnames(exp.data)
  
  ## 构建差异比较矩阵
  contrast.matrix <- makeContrasts(paste0(unique(group_list),collapse = "-"), levels = design)
  
  ## 差异分析
  fit <- lmFit(exp.data, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  tempOutput <- topTable(fit2, coef = 1, n = Inf)
  nrDEG <- na.omit(tempOutput) 
  
  return(nrDEG)
}