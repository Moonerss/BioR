##  TO DO：检测文件目录是否存在以及可写
##
##  Author: jeason zhao
##################################################################

#' @name checkDir
#' @description  检查给定目录是否存在，是否可写
#' @param  inpurDir, 输入目录
#' @returnType 逻辑值向量
#' @return 返回 TRUE/FALSE
#'
checkDir <- function(inpurDirs){
  sapply(inpurDirs, checkDir_core)
}

#' @@ 检查给定目录是否存在，是否可写，注释同上，一次只处理一个目录
checkDir_core <- function(inpurDir){
  if (dir.exists(inpurDir)) {
    # 设置可写
    #writable <- file.access(inpurDir, mode=2)
    # 0 test for existence.
    # 1 test for execute permission.
    # 2 test for write permission.
    # 4 test for read permission.
    Sys.chmod(inpurDir, mode="775", use_umask = FALSE)
  } else {
    cat("create file directory ing......")
    dir.create(inpurDir, showWarnings = TRUE, recursive=TRUE, mode="775")
    Sys.chmod(inpurDir, mode="775", use_umask = FALSE)
  }
}
