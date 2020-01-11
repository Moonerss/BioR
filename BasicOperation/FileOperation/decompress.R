# TODO: 综合解压各种格式的文件
# 
# Author: jeason zhao
###############################################################################

#' @name decompress
#' @description 目前支持zip, tar.gz, tar.bz2, tar, gz
#' @param file 需要进行解压的压缩文件名字
#' @param outdir 解压后的文件路径名
#' @returnType 
#' @return 
decompress <- function(file, outdir=NULL) {
  if(is.null(outdir)){
    outdir <- dirname(file) 
  }
  
  fileType <- basename(file)
  if (grepl("zip$", fileType)) {
    unzip(file, exdir=outdir)
    t1 <- unzip(file, exdir=outdir, list=TRUE)	
    return(as.character(t1$Name))
  } else if(grepl("tar\\.gz$", fileType)) {
    untar(file, exdir=outdir)
    return(untar(file, exdir=outdir, list=TRUE))
  } else if(grepl("tar$", fileType)) {
    untar(file, exdir=outdir)
    return(untar(file, exdir=outdir, list=TRUE))
  } else if(grepl("tar\\.bz2$", fileType)) {
    untar(file, exdir=outdir)
    return(untar(file, exdir=outdir, list=TRUE))
  } else if(grepl("\\.gz$", fileType)) {
    library(R.utils)
    gunzip(file, overwrite=TRUE, remove=FALSE)
    if (outdir != dirname(file)) {
      #解压后的文件地址
      t1 <- gsub("[.]gz$", "", file)
      #拷贝到规定的目录
      file.copy(t1, file.path(outdir, basename(t1)))
      file.remove(t1)
      return(basename(t1))
    } else {
      return("OK")
    }
  }
}


