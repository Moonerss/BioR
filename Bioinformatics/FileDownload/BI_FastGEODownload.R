#' @title GEO下载数据 基于GEOquery
#' @param gseid GSE号
#' @param filepath 文件的存储路径，如："D:/learn"
#' @param log2 是否对下载的表达谱进行log2转化
#' @param origin 是否下载原始cel数据，若origin = TRUE，则使用下载到的原始数据进行表达谱的构建
FastDownloadGEO <- function(gseid, filepath, log2 = F, origin = FALSE) {
  library(GEOquery)
  
  ## 如果id以GSE开头
  if (grepl("GSE",gseid)) {
    # 获取GSE soft信息
    gse <- getGEO(gseid, GSEMatrix = F,destdir = filepath)
    # 判断获取的数据类型(芯片还是测序)
    if (grepl("array",Meta(gse)$type)) {
      # 如果是芯片数据
      cat(paste(gseid," is microarray"))
      # 判断是否下载原始数据
      if (origin == F) {
        # 如果不下载原始数据，则构建文件夹下载表达谱
        dir.name <- file.path(filepath, gseid)
        if(!file.exists(dir.name)) {
          dir.create(dir.name,recursive = T)
        }
        gsets <- getGEO(gseid,getGPL = T, destdir = dir.name)
        # 一个GSE可能存在多个表达谱
        for (i in 1:length(gsets)) {
          gset <- gsets[[i]]
          exprSet <- exprs( gset )  ## “GEOquery”包中的exprs函数用来取出表达矩阵
          pdata <- pData( gset )  ## “GEOquery”包中的pData函数用来取出样本信息
          GPL <- gset@featureData@data  ##探针注释信息
          platform <- gset@annotation ##探针平台
          ## 判断表达谱是否log2标化
          ex <- exprSet
          qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
          LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
          
          ## 判断是否log2转化，并确定是否进行log转化
          if (LogC) {
            if (log2) {
              ex[which(ex <= 0)] <- NaN
              exprSet <- log2(ex)
              print("log2 transform finished")
            }
          } else {
            print("log2 transform not needed")
          }
          ## 导出数据
          temp <- c()
          if (dim(exprSet)[1] > 0) {
            expression.filename <- paste0( gseid, "-" ,platform ,"-expression.MATRIX.txt" )
            write.table(exprSet,file=paste0( filepath,"/", gseid,"/" , expression.filename), sep="\t", quote=F, row.names=T, col.names=T)
            temp <- append(temp,paste0( filepath, "/", gseid,"/" , expression.filename))
          }
          if (dim(pdata)[1] > 0) {
            pData.filename <- paste0( gseid, "-" ,platform ,"-pData.txt" )
            write.table(pdata,file=paste0( filepath, "/",gseid,"/" , pData.filename), sep="\t", quote=F, row.names=F, col.names=T)
            temp <- append(temp,paste0( filepath,"/", gseid,"/" , pData.filename))
          }
          if (dim(GPL)[1] > 0) {
            GPL.filename <- paste0( gseid, "-" ,platform ,"-ProbeAnnotation.txt" )
            write.table(GPL,file=paste0( filepath,"/", gseid,"/" , GPL.filename ), sep="\t", quote=F, row.names=F, col.names=T)
            temp <- append(temp,paste0( filepath,"/", gseid,"/" , GPL.filename))
          }
        }
      } else {
        # 下载原始数据
        getGEOSuppFiles(gseid,makeDirectory = TRUE,baseDir = filepath,fetch_files = TRUE,filter_regex = NULL)
        untar(file.path(file.path(filepath,gseid),paste0(gseid,"_RAW.tar")),exdir=file.path(filepath,gseid))
        library(affy)
        #读入CEL文件
        affy.data <- ReadAffy(celfile.path=file.path(filepath,gseid)) 
        #Background correcting; Normalizing; Calculating Expression
        eset.rma <- rma(affy.data)
        #提取导出表达数据
        exprSet <- exprs(eset.rma) 
        expression.filename <- paste0(gseid,"-raw.expression.MATRIX.txt")
        write.table(exprSet,file=paste0(filepath,"/",gseid,"/",expression.filename),sep="\t",quote=F,row.names=T,col.names=T)
        
        #下载作者处理后的数据提取样本信息以及注释信息
        gsets <- getGEO(gseid,getGPL = T,destdir=file.path(filepath,gseid)) 
        options(stringsAsFactors = F)
        for(i in 1:length(gsets)){
          gset = gsets[[i]]
          exprSet = exprs( gset )  ## “GEOquery”包中的exprs函数用来取出表达矩阵
          pdata = pData( gset )  ## “GEOquery”包中的pData函数用来取出样本信息
          GPL = gset@featureData@data  ##探针注释信息
          platform <- gset@annotation ##探针平台
          temp <- c()
          if(dim(pdata)[1] > 0){
            pData.filename <- paste0( gseid, "-" ,platform ,"-pData.txt" )
            write.table(pdata,file=paste0( filepath,"/", gseid,"/" , pData.filename), sep="\t", quote=F, row.names=F, col.names=T)
            temp <- append(temp,paste0( filepath,"/", gseid,"/" , pData.filename))
          }
          if(dim(GPL)[1] > 0){
            GPL.filename <- paste0( gseid,"-" ,platform , "-ProbeAnnotation.txt" )
            write.table(GPL,file=paste0( filepath,"/", gseid,"/" , GPL.filename ), sep="\t", quote=F, row.names=F, col.names=T)
            temp <- append(temp,paste0( filepath,"/", gseid,"/" , GPL.filename))
          }
        }
      }
    } else if (grepl("sequencing", Meta(gse)$type) & grepl("Expression profiling", Meta(gse)$type)) {
      ## 如果是测序数据，直接下载处理好的数据
      dir.name <- file.path(filepath,gseid)
      if(!file.exists(dir.name)){
        dir.create(dir.name)
      }
      gsets <- getGEO(gseid, getGPL = T, destdir = dir.name)
      for (i in 1:length(gsets)) {
        gset <-  gsets[[i]]
        exprSet <- exprs(gset) ## “GEOquery”包中的exprs函数用来取出series matrix矩阵
        if( dim(exprSet)[1] > 0 ){
          expression.filename <- paste0(gseid, "-expression.MATRIX.txt")
          write.table(exprSet,file=paste0( filepath,"/", gseid,"/" , expression.filename), sep="\t", quote=F, row.names=T, col.names=T)
          cat(paste0(gseid," processed RNA-seq raw data downloaded"))
        }else{
          cat(paste0(gseid," only RNA-seq raw data，please download according to aspera or SRAtool-kit")) ## RNA-seq原始数据就不下载了
        }
      }
    }
  }
}

#' @usage FastDownloadGEO("GSE31519","D:",origin = T) FastDownloadGEO("GSE77445","D:",origin = F)



#---------------------- 判断芯片表达谱是否需要log2标化 ---------------------#
#ex <- exprSet
#qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
#LogC <- (qx[5] > 100) ||
#  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
#  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

#if (LogC) { 
#  ex[which(ex <= 0)] <- NaN
#  exprSet <- log2(ex)
#  print("log2 transform finished")
#}else{
#  print("log2 transform not needed")
#}




