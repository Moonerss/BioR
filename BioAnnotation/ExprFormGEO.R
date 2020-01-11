######################################
# TO DO：完成GEO microarray的注释
# 
# Author：jeason zhao
######################################

#' @name annoExprs.fromGEOMatrix
#' @description probe matrix convert to gene expression matrix
#' @param t.matrix 必须是matrix类型，行名是探针ID，
#' @param GPL.matrix 必须是matrix，对应proID和基因ID必须是charater，顺序是proID、基因ID
#' @param method 一个字符串：c('mean','median')对应多个探针指到同一个基因上，对应基因表达值使用均值还是中值
#' @returnType matrix
#' @return 返回基因表达矩阵

annoExprs.fromGEOMatrix <- function(t.matrix, GPL.matrix, method="median")
{
  # GPL.matrix
  # t.id <- "ID"
  # t.geneid <- "SPOT_ID"
  # t.anno <- IPD(GPL)
  t.anno <- GPL.matrix
  
  
  
  #映射到多个基因或者映射不明确，直接删掉
  t.index <- rep(TRUE, nrow(t.anno))
  if(any(grepl(" ", as.character(t.anno[,2])))){
    #有些geneid包含多个基因，默认认为是空格隔开，如果有则删掉
    t.index <- sapply(as.character(t.anno[,2]), function(x) length(unique(unlist(strsplit(x, " ")))))==1;
  }
  #不为空
  t.index <- t.index & (t.anno[,2] != "")
  #不为NA
  t.index <- t.index & (!is.na(t.anno[,2]))
  
  t.anno <- t.anno[t.index, ]
  
  #应该全部都能match上
  temp1 <- t.anno[match(rownames(t.matrix), t.anno[, 1]), 2];
  
  rownames(t.matrix) <- temp1;
  t.result <- meanMultipleSameGene(t.matrix, method=method);
  
  if(all(is.na(rownames(t.matrix)))){
    stop("所有基因名字均为NA，不同情况导致GPL转换失败:\n
         1)完全无法注释的GPL，直接删除该数据；\n
         2)检查该数据探针名是否与GPL中的一致，存在大小写问题")
  }
  
  return(t.result)
}

#' @name meanMultipleSameGene
#' @description  多行对应于一个基因，直接取他们的均值或者中值
#' @param t.matrix 必须是matrix类型
#' @param method 均值还是中值，暂时只支持均值
#' @returnType 
#' @return 
meanMultipleSameGene <- function(t.matrix, method="median") {
  
  #考虑到计算时间很慢，因此，把矩阵分成两部分
  #得到唯一性的名字列表
  t.name <- sort(unique(rownames(t.matrix)))
  
  #取出出现至少两次的names
  t.two.name <- names(which(table(rownames(t.matrix)) > 1));
  #如果没有出现两次的名字，直接返回当前表达矩阵
  if (length(t.two.name) == 0) {
    return(t.matrix)
  }
  t.two.matrix <- t.matrix[rownames(t.matrix) %in% t.two.name,]
  t.one.matrix <- t.matrix[!(rownames(t.matrix) %in% t.two.name),]
  
  if (method=="mean") {
    t1 <- by(t.two.matrix, rownames(t.two.matrix), colMeans, na.rm = TRUE)		
  }
  if (method=="median") {
    t1 <- by(t.two.matrix, rownames(t.two.matrix), function(x) apply(x, 2, median, na.rm = TRUE))		
  }	
  t.result <- do.call(rbind, t1)
  
  t.result <- rbind(t.result, t.one.matrix);
  
  t.result <- t.result[t.name,]
  return(t.result)
}