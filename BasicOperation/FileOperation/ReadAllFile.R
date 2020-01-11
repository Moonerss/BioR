#########################################
# TODO: 批量读取一个文件夹下的文件
#
# Author: Jeason Zhao
#########################################

#########################################
##MAIN:批量读取文件夹下的txt文件
##Ref:文件必须是同一种文件类型
#########################################

#' 读取同一目录下的所有txt文件
#' @param path 文件所在的路径

readAlltxt <- function(path) {
  ##获取该路径下的文件名
  fileNames <- dir(path)
  filePath <- sapply(fileNames, function(x) { 
                     paste(path,x,sep='/')
  })   ##生成读取文件路径
  data <- lapply(filePath, function(x) {
    read.table(x, header=T, sep = "\t", stringsAsFactors = FALSE, quote = "\"")
  })  ##结果生成一个list
}




#' 读取同一目录下的所有txt文件并合并成一个文件
#' 这里一般文件的内容格式是相似的，如：不同样本的突变数据
#' @param path 文件所在的路径

readAlltxt.Merge <- function(path) {
  fileNames <- dir(path);  ##获取该路径下的文件名
  filePath <- sapply(fileNames, function(x) { 
    paste(path,x,sep='/')
  }) ##生成读取文件路径
  merge.data <- read.table(file = filePath[1], header = T, sep = "\t", stringsAsFactors = FALSE, quote = "\"");
  for(i in 2:length(filePath)) {
    new.data <- read.table(file = filePath[i], header = T, sep = "\t", stringsAsFactors = FALSE, quote = "\"");
    merge.data <- rbind(merge.data, new.data)
  }
  return(merge.data)
}