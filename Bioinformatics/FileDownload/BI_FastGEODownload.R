#' @title GEO下载数据 基于GEOquery
#' @param GSEid
#' @param filepath eg. "D:/learn"
#' @param origin 是否下载原始cel数据

FastDownloadGEO <- function(GSEid,
                        filepath,
                        origin=F){
  library(GEOquery)
  
  #如果是GSE开头
  if(grepl("GSE",GSEid)){  
    #获取GSE soft信息
    gse <- getGEO(GSEid,GSEMatrix = F,destdir=filepath)
    # 判断是芯片数据还是测序数据 
    if(Meta(gse)$type == "Expression profiling by array" || Meta(gse)$type == "Non-coding RNA profiling by array" || grepl("array",Meta(gse)$type))
      { 
      # 如果是芯片数据并且存在
      cat(paste(GSEid," is microarray"))
        ## 芯片非原始数据下载
        if(origin == F){
          dir.name <- file.path(filepath,GSEid)
          if( !file.exists(dir.name) ){
            dir.create(dir.name,recursive = T)
          }
          gsets <- getGEO(GSEid,getGPL = T,destdir=dir.name) #下载数据 
          options(stringsAsFactors = F)
		  ## 一个GSE可能存在多个表达矩阵
		  for(i in 1:length(gsets)){
			gset <- gsets[[i]]
			exprSet <- exprs( gset )  ## “GEOquery”包中的exprs函数用来取出表达矩阵
			pdata <- pData( gset )  ## “GEOquery”包中的pData函数用来取出样本信息
			GPL <- gset@featureData@data  ##探针注释信息
			platform <- gset@annotation ##探针平台
			## 判断表达谱是否log2标化
			ex <- exprSet
			qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
			LogC <- (qx[5] > 100) ||
			(qx[6]-qx[1] > 50 && qx[2] > 0) ||
			(qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
			if(LogC){ 
				ex[which(ex <= 0)] <- NaN
				exprSet <- log2(ex)
				print("log2 transform finished")
			}else{
				print("log2 transform not needed")
			} 
			temp <- c()
			if(dim(exprSet)[1] > 0){
				expression.filename <- paste0( GSEid, "-" ,platform ,"-expression.MATRIX.txt" )
				write.table(exprSet,file=paste0( filepath,"/", GSEid,"/" , expression.filename), sep="\t", quote=F, row.names=T, col.names=T)
				temp <- append(temp,paste0( filepath, "/", GSEid,"/" , expression.filename))
			}
			if(dim(pdata)[1] > 0){
				pData.filename <- paste0( GSEid, "-" ,platform ,"-pData.txt" )
				write.table(pdata,file=paste0( filepath, "/",GSEid,"/" , pData.filename), sep="\t", quote=F, row.names=F, col.names=T)
				temp <- append(temp,paste0( filepath,"/", GSEid,"/" , pData.filename))
			}
			if(dim(GPL)[1] > 0){
				GPL.filename <- paste0( GSEid, "-" ,platform ,"-ProbeAnnotation.txt" )
				write.table(GPL,file=paste0( filepath,"/", GSEid,"/" , GPL.filename ), sep="\t", quote=F, row.names=F, col.names=T)
				temp <- append(temp,paste0( filepath,"/", GSEid,"/" , GPL.filename))
			}
		  }
        }else{
          #下载Cel原始文件
          getGEOSuppFiles(GSEid,makeDirectory = TRUE,baseDir = filepath,fetch_files = TRUE,filter_regex = NULL)
          untar(file.path(file.path(filepath,GSEid),paste0(GSEid,"_RAW.tar")),exdir=file.path(filepath,GSEid))
          library(affy)
          #读入CEL文件
          affy.data <- ReadAffy(celfile.path=file.path(filepath,GSEid)) 
          #Background correcting; Normalizing; Calculating Expression
          eset.rma <- rma(affy.data)
          #提取表达
          exprSet <- exprs(eset.rma) 
          expression.filename <- paste0(GSEid,"-raw.expression.MATRIX.txt")
          write.table(exprSet,file=paste0(filepath,"/",GSEid,"/",expression.filename),sep="\t",quote=F,row.names=T,col.names=T)
          
          #下载作者处理后的数据提取样本信息以及注释信息
          gsets <- getGEO(GSEid,getGPL = T,destdir=file.path(filepath,GSEid)) 
          options(stringsAsFactors = F)
		  for(i in 1:length(gsets)){
			gset = gsets[[i]]
			exprSet = exprs( gset )  ## “GEOquery”包中的exprs函数用来取出表达矩阵
			pdata = pData( gset )  ## “GEOquery”包中的pData函数用来取出样本信息
			GPL = gset@featureData@data  ##探针注释信息
			platform <- gset@annotation ##探针平台
			temp <- c()
			if(dim(pdata)[1] > 0){
				pData.filename <- paste0( GSEid, "-" ,platform ,"-pData.txt" )
				write.table(pdata,file=paste0( filepath,"/", GSEid,"/" , pData.filename), sep="\t", quote=F, row.names=F, col.names=T)
				temp <- append(temp,paste0( filepath,"/", GSEid,"/" , pData.filename))
            }
            if(dim(GPL)[1] > 0){
				GPL.filename <- paste0( GSEid,"-" ,platform , "-ProbeAnnotation.txt" )
				write.table(GPL,file=paste0( filepath,"/", GSEid,"/" , GPL.filename ), sep="\t", quote=F, row.names=F, col.names=T)
				temp <- append(temp,paste0( filepath,"/", GSEid,"/" , GPL.filename))
			}
		  }
          
        }
      }
    }else if(grepl("sequencing",Meta(gse)$type) & grepl("Expression profiling",Meta(gse)$type) ){ ##是否是测序数据
      ## 首先判断是否提供处理好的数据
      dir.name <- file.path(filepath,GSEid)
      if(!file.exists(dir.name)){
        dir.create(dir.name)
      }
      gsets <- getGEO(GSEid,getGPL = T,destdir=dir.name) 
      options(stringsAsFactors = F)
	  for(i in 1:length(gsets)){
		gset = gsets[[i]]
        exprSet = exprs(gset) ## “GEOquery”包中的exprs函数用来取出series matrix矩阵
        if( dim(exprSet)[1] > 0 ){
			expression.filename <- paste0(GSEid, "-expression.MATRIX.txt")
			write.table(exprSet,file=paste0( filepath,"/", GSEid,"/" , expression.filename), sep="\t", quote=F, row.names=T, col.names=T)
			cat(paste0(GSEid," processed RNA-seq raw data downloaded"))
        }else{
			cat(paste0(GSEid," only RNA-seq raw data，please download according to aspera or SRAtool-kit")) ## RNA-seq原始数据就不下载了
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




