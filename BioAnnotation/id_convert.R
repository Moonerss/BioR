#' @TODO 基因ID转换
#' @param input 输入需要转换的基因，字符型向量
# > head(input)
# [1] "RRM2"   "HMMR"   "NUSAP1" "UBE2C"  "TOP2A"  "AURKA" 
#' @param from 输入的类型
#' @param to 转换的类型
#' @param database 使用的注释数据库
#' @returnType data.frame
#' @return 转换后的矩阵，两列，第一列为from，第二列为to
# > res
#      SYMBOL ENTREZID
# 1      RRM2     6241
# 2      HMMR     3161
# 3    NUSAP1    51203
# 4     UBE2C    11065
# 5     TOP2A     7153
#
#' @author ZGX
#
id_convert <- function(input=NULL, from="SYMBOL", to="ENTREZID", database="org.Hs.eg.db"){
  library(org.Hs.eg.db)
  library(clusterProfiler)
  res <- bitr(input, fromType=from, toType=to, OrgDb=database)
  # fromType和toType的取值范围如下
  #   > columns(org.Hs.eg.db)
  #  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
  #  [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
  # [11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"        
  # [16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
  # [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"     
  # [26] "UNIPROT" 
  return(res)
  # 注意，也可以使用 bioMart 对基因进行ID转换，见 https://www.jianshu.com/p/3a0e1e3e41d0
}
