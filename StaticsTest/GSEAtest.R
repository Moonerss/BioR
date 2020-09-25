#' @@ GSEAtest
#' @param  gene.rank 给定一个排秩向量，向量值代表排序的依据比如差异表达fold change，名字为基因
#' @param  signatureSet.df 待富集的signature 基因集合，dataframe格式，第一列为signature名称，第二列为基因
#' head(signatureSet.df)
#   SetAnno                ENSEMBL
# 1 leukocyte_infiltration ENSG00000123384
# 2 leukocyte_infiltration ENSG00000115457
# 3 leukocyte_infiltration ENSG00000105855
#' @param  minGSSize 映射到排秩列表的基因数目不少于，默认为10
#' @param  nPerm 随机次数，默认为1000
#' @param  pAdjustMethod p值矫正方法，默认为"BH"
#' @param  pvalueCutoff p的阈值，如果设为1，则输出每个signature的富集结果
#' @param  plotSignature 可以自己指定画哪个signature，默认等于all，所有的都画出
#' @param  output.file 绘图结果存储文件名
#' @returnType 
#' @return 
#' 
#' @author 
GSEAtest <- function(gene.rank, signatureSet.df, 
                     nPerm = 1000, minGSSize = 10, pvalueCutoff = 1, pAdjustMethod = "BH",
                     plotSignature = "all", output.file){
  
  library(clusterProfiler)
  library(enrichplot)
  
  ###1、GSEA富集
  GSEAres <- GSEA(r, nPerm = nPerm, minGSSize = minGSSize, pvalueCutoff = pvalueCutoff, 
                  pAdjustMethod = pAdjustMethod, TERM2GENE = signatureSet.df)
  
  ###2、画图
  pdf(file = output.file)
  
  if(plotSignature == "all"){
    for(sig in unique(GSEAres$ID)) {
      p <- gseaplot2(GSEAres, geneSetID = which(GSEAres$ID == sig), title = plotSignature, pvalue_table = T)
      print(p)
    }
  }else{
    p <- gseaplot2(GSEAres, geneSetID = which(GSEAres$ID == plotSignature), title = plotSignature, pvalue_table = T)
    print(p)
  }
  dev.off()
  
  return(GSEAres)
}


