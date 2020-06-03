####################################################################################
# 整体功能简介：
# 使用tximport包将kallisto/salmon定量的转录本表达转换为基因表达水平
###################################################################################

###################################################################################
# transcript2gene.quant
####原理
#' tximport package支持salmon和kallisto进行转录本和基因表达水平的转换
#' tximport只是导入kallisto/salmon定量的转录本TPM和count值，将同一个基因转录本异构体加和，计算出基因水平的TPM和count

####支持证据
#' salmon官网上有直接说明可以使用tximport包将转录本水平表达转换为基因水平表达
#' kallisto并没有直接的说明可以直接用tximport进行转录本和基因水平表达的转换，kallisto的开发人员推荐使用他们开发的专门适用kallisto下游分析sleuth软件进行处理
#' 但是确实有许多发表的文章使用了kallisto+tximport衔接运行

#' @param quantDir kallisto/salmon定量结果目录，如"/home/longzl2017/dataScience/PipeLine/RNA-seq/test/fastq/kresult"， "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/fastq/sresult"
#' @param outDir 结果输出目录，如"/pub6/temp/immunophenoscore/Data/28552987" 
#' @param type 使用软件类型。目前只支持salmon和kallisto 
#' @param is.RefGencode TRUE/FALSE，kallisto/salmon建立index时使用的参考转录组是否是gencode，推荐使用gencode转录组。如果为真，则tx2gene可以为NULL。但如果用户给定了tx2gene，即使是gencode注释类型，程序也会直接采用用户给定的tx2gene
#' @param tx2gene data.frame，转录本和基因的对应关系，要求第一列为转录本id，第二列为基因id，要求和参考转录组使用的版本一致，否则可能会出现转录本和基因的版本不一致情况，导致无法匹配。
#' 如果gencodeAnnotation为TRUE时，则可以为NULL
#' @param countsFromAbundance 设置从abundance估计count的方式。支持三种方式"no"、"scaledTPM"和"lengthScaledTPM"，通常情况下选择"no"。
#' tximport默认为no，即直接从定量软件中加和同一个基因转录本count值；
#' scaledTPM，基于library size和定量软件计算出TPM重新估计count；
#' lengthScaledTPM，基于library size、所有样本中平均转录本长度和定量软件计算出TPM重新估计count

#' @returnType list
#' @return gene.expression.matrix 包含3个matrix对象：TPM, counts, length。
#' TPM为基因表达TPM矩阵，counts为基因表达count矩阵，length为每个基因的平均转录本长度

#' @note 
#' 当前使用tximport package版本为1.6.0
#' scaledTPM、lengthScaledTPM两种方式，基于TPM重新估计每个基因的count，而不是直接加和软件估计的转录本count值，这种方式能够帮助矫正文库差异和转录本长度

#' example
# 非gencode注释 --- kallisto
#' load(file = "/pub6/temp/xlw/kallisto_linux-v0.44.0/Homo_sapiens.ENSG2ENST2.RData")
#' tx2gene <- data.frame(transcript_id = Homo_sapiens.ENSG2ENST$transcript_id, gene_id = Homo_sapiens.ENSG2ENST$Ensembl_id)，格式可以参考 https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
#' transcript2gene.quant(quantDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/fastq/kresult", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/fastq/kresult", type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "no")

# gencode注释 --- kallisto
#' transcript2gene.quant(quantDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq/kresult", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq/kresult", type = "kallisto", is.RefGencode = T, tx2gene = NULL, countsFromAbundance = "no")

# gencode注释 --- salmon
#' transcript2gene.quant(quantDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq/sresult", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq/sresult", type = "salmon", is.RefGencode = T, tx2gene = NULL, countsFromAbundance = "no")


transcript2gene.quant <- function(quantDir, outDir, type = c("salmon", "kallisto"), is.RefGencode = FALSE, tx2gene, countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM")){
  
  if(!dir.exists(quantDir)){
    stop(quantDir, "not exist!please check...")
  }
  
  # 判断type参数输入是否正确
  if(!(type %in% c("salmon", "kallisto"))){
    stop("type参数输入格式不符!")
  }
  
  if(type == "kallisto"){
    pattern <- "abundance.tsv"
  }
  if(type == "salmon"){
    pattern <- "quant.sf"
  }
  
  quantFiles <- list.files(quantDir, pattern=paste0("^", pattern, "$"), full.names=TRUE, recursive=T)
  if(length(quantFiles) == 0){
    stop(quantDir, "没有匹配相应类型的文件，请确认目录下是否有文件或输入的type参数和定量文件软件来源是否一致!")
  }
  
  sampleName <- list.files(quantDir, pattern=paste0("^", pattern, "$"), full.names=F, recursive=T)
  sampleName <- gsub(paste0("/", pattern), "", sampleName)
  names(quantFiles) <- sampleName
  
  #由于使用gencode中定量出quant.sf或abundance.tsv文件包含转录本和基因名字信息，且样本间是一致的，所以可以利用第一列信息来构建tx2gene
  if(is.RefGencode){
    # 如果用户给定了tx2gene，即使为gencode注释类型也采用用户给定的tx2gene
    if(is.null(tx2gene)){
      geneInfo <- read.table(file = quantFiles[1], sep = "\t", stringsAsFactors = F, header = T)
      geneInfo <- geneInfo[,1]
      infoMatrix <- sapply(geneInfo, function(x){
        attr <- unlist(strsplit(x, "\\|"))
        return(attr)
      })
      infoMatrix <- t(infoMatrix)
      tx2gene <- data.frame(quant.id = rownames(infoMatrix), gene.id = infoMatrix[, 2])  
    }
  }else{
    if(is.null(tx2gene)){
      stop("如果不是gencode注释体系，需要给定tx2gene文件")
    }
  }
  
  library(tximport)
  
  # 判断countsFromAbundance参数输入是否正确
  if(!(countsFromAbundance %in% c("no", "scaledTPM", "lengthScaledTPM"))){
    stop("countsFromAbundance参数输入格式不符!")
  }
  
  # txOut = FALSE表示估计基因水平的表达
  gene.expression.matrix <- tximport(files = quantFiles, type = type, tx2gene = tx2gene, txOut = FALSE, countsFromAbundance = countsFromAbundance)
  gene.TPM.matrix <- gene.expression.matrix$abundance
  gene.count.matrix <- gene.expression.matrix$counts
  gene.length.matrix <- gene.expression.matrix$length
  
  gene.expression.matrix <- list(TPM = gene.TPM.matrix, count = gene.count.matrix, length = gene.length.matrix)
  save(gene.expression.matrix, file = file.path(outDir, "gene.expression.matrix.RData"))
}


