#' @@ DESeq2Diff
#' 基于DESeq2包计算RNA-seq差异表达
#' @param  countData 输入的原始read counts矩阵。行为基因，列为样本。
#' @param  colData 输入的样本注释矩阵。colData的行名与countData的列名一致。
# > head(colData)
#          condition  condition2
# control1   control  ttt
# control2   control  sss
# treat1       treat  sss
# treat2       treat  sss
#' @param  classCol 表示将colData中的哪一列作为样本分类标签，给出列名
#' @param  contrast 设定在计算差异表达的时候如何进行比较
# 当classCol只有一个条件时，contrast是长度为3的字符型向量，分别指定
# the name of a factor in the design formula, 
# the name of the numerator level for the fold change,
# and the name of the denominator level for the fold change (simplest case)
# 例如：
# contrast=c("condition", "A", "B") 结果中的foldchange = A/B
# contrast=c("condition", "B", "A") 结果中的foldchange = B/A
#' @param pval_cutoff padj阈值
#' @param pAdjustMethod p值矫正方法，"BH", "fdr", "none"
#' @param log2fc_cutoff log2-transformed fold change阈值
#' @param return_df 是否以data.frame返回
#' @returnType 
#' @return 
#' 
#' @author 

DESeq2Diff <- function(countData=NULL, colData=NULL, classCol="condition", contrast=c("condition", "B", "A"), 
                       pval_cutoff=0.05, pAdjustMethod=c("BH", "fdr", "none"), log2fc_cutoff=1, return_df=TRUE){
  library(DESeq2)
  # 构建DESeq2对象所需的design formula
  design_formula <- as.formula(paste("~", classCol))
  # 构建DESeq2对象
  dds_obj <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = design_formula)
  # 基于DESeq2计算差异
  dds <- DESeq(dds_obj)
  # 提取结果
  res <- results(dds, contrast = contrast, pAdjustMethod = pAdjustMethod)
 
  # 过滤 padj
  if(!is.null(pval_cutoff)){
    res <- subset(res, padj < pval_cutoff)
  }
  # 过滤log2-transformed fold change
  if(!is.null(log2fc_cutoff)){
    res <- subset(res, abs(log2FoldChange) > log2fc_cutoff)
  }
  #排序结果矩阵
  if(dim(res)[1]>0){
    res <- res[order(res$padj), ]
    # 是否以data.frame返回
    if(return_df){
      return(as.data.frame(res))
    }else{
      return(res)
    }
  }else{
    return(NA)
  }
}

