###################################################
# TO DO：FPKM、RPKM、Count matrix转化为TPM
# 
# Author: jeason zhao
###################################################


#' @description 输入基因水平或转录本水平的表达矩阵（FPKM、RPKM、Count），返回转化成TPM之后的表达矩阵
#' @param  Exp 表达矩阵，matrix，行为基因列为样本
#' @param  matrices Exp的表达测度，字符，默认是Fpkm，可选参数为count、Fpkm、Rpkm
#' @param  log Exp是否经过log转化，逻辑值，默认为F
#' @param  logbase 如果Exp经过了log转化，数值，log底是多少，默认为NULL
#' @param  GeneExonLen 基因外显子长度信息，字符向量，并且具有names()信息,默认为空
#' @returnType matrix
#' @return TPMExp，且没有经过log转化

ConvertToTPM <- function(Exp, matrices = "Fpkm", GeneExonLen = NULL, log = F, logbase = NULL){
  ## 确保基因表达谱中的基因和基因长度文件中的基因一致，否则取交集
  if (!is.null(GeneExonLen)) {
    Inconsistent <- setdiff(rownames(Exp), names(GeneExonLen)) #表达矩阵和基因长度信息中的不一致的基因  
    if(length(Inconsistent) > 0){
      comGene <- intersect(rownames(Exp), names(GeneExonLen))
      Exp <- Exp[comGene, ]
      GeneExonLen <- GeneExonLen[comGene]
    }
    if(length(Inconsistent) == 0){
      Exp <- Exp #如果基因表达谱中的基因和基因长度文件中的基因一致，则不做任何处理
      GeneExonLen <- GeneExonLen[rownames(Exp)]
    }
  } else {
    stop("GeneExonLen is NULL!")
  }
  
  #确保Exp是没有经过log转化
  if(log==T){
    Exp <- logbase^Exp -1 #如果表达矩阵经过了log转化，则逆对数转回去
  }
  
  #caculate gene length
  if(is.null(names(GeneExonLen))){
    GeneExonLen <- GeneExonLen[rownames(Exp)] #如果基因长度没有按照表达谱中基因的顺序排列，则变换顺序
  }
  
  #convert count to TPM
  if(matrices == "count"){
    TPMExp <- apply(Exp, 2, function(x){
      CountToTpm(counts = x, lengths = GeneExonLen)
    })
  }
  
  #convert Fpkm to TPM
  if(matrices == "Fpkm"){
    TPMExp <- apply(Exp, 2, function(x){
      FpkmOrRpkmToTpm(FpkmOrRpkm = x)
    })
  }
  
  #convert Rpkm to TPM
  if(matrices == "Rpkm"){
    TPMExp <- apply(Exp, 2, function(x){
      FpkmOrRpkmToTpm(FpkmOrRpkm = x)
    })
  }
  
  return(TPMExp)
  
}

#---convert count to TPM----
CountToTpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

#---convert Fpkm/Rpkm to TPM---
FpkmOrRpkmToTpm <- function(FpkmOrRpkm) {
  FpkmOrRpkm*1e6/sum(FpkmOrRpkm) 
  #TPM = ( FPKM / sum of FPKM over all genes/transcripts) * 10^6
  #TPM = ( RPKM / sum of RPKM over all genes/transcripts) * 10^6
}

#' @description  输入从biomart中下载的基因转录本的信息，返回每个基因的外显子长度
#' @param  transcriptInfor  从biomart中下载的基因转录本的信息，其中第一列为基因名，第二列为转录本名，第三列为外显子起始位置，第四列为外显子终止位置
#' @return TranLen，每个基因的外显子长度
#' 
#----获得外显子长度-----
getExonLength <- function(transcriptInfor){
  #输入从biomart中下载的基因转录本的信息，其中第一列为基因名，第二列为转录本名，第三列为外显子起始位置，第四列为外显子终止位置
  #输出每个基因的外显子长度
  ExonLen <- transcriptInfor[,"Exon.region.end..bp."]- transcriptInfor[,"Exon.region.start..bp."] #外显子长度
  TranLen <- tapply(ExonLen, transcriptInfor[,"Gene.stable.ID"], sum) #每个基因所有外显子的长度
  return(TranLen)
}

