##  TO DO: 使用R进行并行运算时需要的内存计算操作
##
##  Author: jeason zhao
################################################################


#' @name sysmem
#' @description 返回系统总内存数
#' @returnType 
#' @return 
#' 
sysmem <- function() {
  t1 <- system("cat /proc/meminfo", intern=TRUE);
  t1 <- gsub("[[:alpha:][:punct:][:space:]]", "", t1[1])
  return(as.numeric(t1))
}


#' @name memR
#' @description 当前可用内存为多少
#' @returnType 
#' @return 
#' 
memR <- function() {
  bit <- 8L * .Machine$sizeof.pointer
  if (!(bit == 32L || bit == 64L)) {
    stop("Unknown architecture", call. = FALSE)
  }
  
  node_size <- if (bit == 32L) 28L else 56L
  
  usage <- gc()
  sum(usage[, 1] * c(node_size, 8)) / (1024 ^ 2)
}


#' @name memLinux
#' @description 返回当前R在linux中占用的内存， Mb
#' @return 
memLinux <- function(){
  library(fork)
  temp <- as.numeric(system(paste("ps -p ", getpid(), " -o rss", sep=""), intern=TRUE)[2])
  return(temp/1024)
}


#' @name autoCore
#' @description 自动检测适合的核用于并行运算
#' @param full : TRUE表示使用所有核
#' @returnType 
#' @return 
autoCore <- function(full=FALSE) {
  library(parallel)
  if (full) {
    mc.cores <- detectCores()
    return(mc.cores)
  }
}


