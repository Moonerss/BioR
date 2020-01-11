##  TO DO: 将list转换成其他格式
##
##  Author：jeason zhao
#########################################

#' @name list2matrix
#' @description 将list转换成0,1矩阵，判断每个元素是否在list中出现
#' @param t.list 
#' @return 返回0,1矩阵
list2matrix <- function(t.list){
  t1 <- unique(unlist(t.list));
  t2 <- matrix(0, ncol=length(t.list), nrow=length(t1));
  rownames(t2) <- t1;
  colnames(t2) <- names(t.list);
  
  for(i in 1:length(t.list)){
    t2[,i] <- as.numeric(rownames(t2) %in% t.list[[i]]);
  }
  return(t2)
}