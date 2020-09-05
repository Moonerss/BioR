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
#' @param pattern 正则表达式，用于匹配特定文件
#' @param merge 是否合并读取的文件，默认为TRUE，合并文件；否则返回一个list
#' @param verbose 输出额外提示信息
#' @param ... fread函数的参数

read_all_file <- function(path, pattern = NULL, merge = TRUE, verbose = TRUE,...) {
    if (!requireNamespace("data.table")) {
        stop("Your need install the `data.table` package")
    }
    if (verbose) message("=> Starting")
    filenames <- dir(path, pattern = pattern)
    filepath <- sapply(filenames, function(x) {paste(path,x,sep='/')})
    if (length(filepath) > 1) {
        if (merge) {
            if (verbose) message("==> Reading ", filenames[1])
            merge_data <- data.table::fread(file = filepath[1], header = T, sep = "\t", stringsAsFactors = FALSE, ...)
            for (i in 2:length(filepath)) {
                if (verbose) message("==> Reading ", filenames[i])
                new_data <- data.table::fread(file = filepath[i], header = T, sep = "\t", stringsAsFactors = FALSE, ...);
                merge_data <- rbind(merge_data, new_data)
            }
            res <- merge_data
        } else {
            if (verbose) message("==> Reading files ...")
            res <- data <- lapply(filepath, data.table::fread, ...)
        }  
    } else {
        if (verbose) message("==> Reading file ...")
        res <- data.table::fread(file = filepath, ...)
        res <- as.data.frame(res)
    }
    
    if (verbose) message("=> Done!")
    return(res)
}