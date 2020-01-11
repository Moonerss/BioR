##  TO DO: SRA文件转fastq
##
##  Author: jeason zhao
####################################


###################################################################################
#函数名称：convertSRA2FASTQ
#' @param filePath: 待转换的.sra文件存放路径及文件名（仅用于每次处理单个SRA文件转换时）如：/IData/SingleCell/scRNAseq/Brain/SRP092584/SRR4919521/SRR4919521.sra
#' @param outDir: 数据分析结果存储目录（切记输出目录定义的时候，正确：/pub5/xiaoyun/OurData，错误写法：/pub5/xiaoyun/OurData/），会自动建立fastq文件夹目录用于存放fastq文件，如outDir = "/pub5/xiaoyun/OurData"，则会在outData下建立fastq文件夹
# 在outDir为NULL情况下，对于批量sra文件和单个sra文件处理，会有不同的处理：
# 多sra文件处理，会自动在dirPath路径下建立名为fastq文件夹存储转换后的结果文件
# 单sra文件处理，会自动在sra文件所在目录下生成转换的fastq文件
#' @param dirPath: 待转换的.sra文件集的文件路径（仅用于SRA文件批量转换时）
#' @param pattern: 模式匹配sra文件，默认为".sra$"
#' @param pairEND: TURE/FALSE，TURE则表示为双末端测序，FALSE则表示为单末端测序（单个sra文件可能存储了pair-end的两个配对的read，因此，该参数可以从SRA文件中分别提取出两个_1与_2的FASTQ文件）
#' @param threads: 并行线程数目。一般默认为20，适合批量sra文件转换
#' @return 返回转换后的fastq.gz的文件

#' @Note:
# 函数 convertSRA2FASTQ 使用的内核工具SRAtoolkit为最新版本 2.8.2
# 函数内部内嵌了glibc/python-3.6.5/parallel-fastq-dump-master软件的路径，都在/pub5/xiaoyun/BioSoftware 目录下，运行时注意检查是否存在
# 批量处理sra文件时，同时会在outDir目录下产生一个sra2fastq_logFile.txt文件，对于每个样本通过检查是否存在转换的fastq文件来判断是否转换成功。

#### convertSRA2FASTQ实现了并行fastq-dump
# 基础的NCBI fastq-dump并不支持并行，大文件处理时运行速度极慢，即使你有多余的CPU，也无法加快速度。
# 由于fastq-dump有选项（-N和-X）来查询sra文件的特定范围，parallel-fastq-dump(https://github.com/rvalieris/parallel-fastq-dump)工具通过将一个sra切割为多份，分配给多个线程，并行运行多个fastq-dump，最终将结果合并
# 经过本人测试，parallel-fastq-dump和NCBI fastq-dump得到的结果文件大小一致，转换的spots数目也是相同的

##example
# single sra
# convertSRA2FASTQ(filePath = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single_sra/SRR1821348.sra", pairEND = F, threads = 2)
# multi-sra paired
# convertSRA2FASTQ(outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired", dirPath = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired", pairEND = T, threads = 10)
# multi-sra single
# convertSRA2FASTQ(outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single", dirPath = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single", pairEND = F, threads = 10)


convertSRA2FASTQ <- function(filePath = NULL, outDir = NULL, dirPath = NULL, pattern = ".sra$", threads = 20, pairEND = FALSE){
  
  # 导入相关的环境变量
  Sys.setenv(LD_LIBRARY_PATH = "/pub5/xiaoyun/BioSoftware/glibc/lib:$LD_LIBRARY_PATH")
  Sys.setenv(PATH = "/pub5/xiaoyun/BioSoftware/Python-3.6.5/bin:/pub5/xiaoyun/BioSoftware/sratoolkit-2.8.2/bin:$PATH")
  Sys.setenv(PYTHONPATH = "/pub5/xiaoyun/BioSoftware/Python-3.6.5/lib/python3.6:$PYTHONPATH")
  
  # 运行命令
  SRA.dump.fastq <- "/pub5/xiaoyun/BioSoftware/parallel-fastq-dump-master/parallel-fastq-dump"
  
  # 批量处理sra文件
  if(!is.null(dirPath)){
    # 匹配所有的sra文件
    t.files <- list.files(dirPath, pattern = pattern, full.names = TRUE, recursive = TRUE)
    
    # 确定样本名字，设定.sra前面的字段为样本名字 
    sampleNames <- list.files(dirPath, pattern = pattern, full.names = TRUE, recursive = TRUE)
    sampleNames <- gsub(".*/", "", sampleNames)
    sampleNames <- gsub(".sra$", "", sampleNames)
    
    # 创建目录
    if(is.null(outDir)){
      outDir <- dirPath
    }
    
    # 建立日志文件，便于后续处理转换有问题的样本
    sra2fastq_log_file <- file.path(outDir, "sra2fastq_logFile.txt")
    
    outDir <- file.path(outDir, "fastq")
    if(!dir.exists(outDir)){
      dir.create(outDir)
    }
    
    logout <- c()
    for(i in 1:length(t.files)){
      if(!pairEND){
        t.command <- paste(SRA.dump.fastq, "-s", t.files[i], "-t", threads, "--gzip -O", outDir, sep = " ")
        system(t.command)
        fastqFile <- paste0(sampleNames[i], ".fastq.gz")
        fastqFile <- file.path(outDir, fastqFile)
        
        # 通过检查转换的fastq文件是否存在，来判断转换是否成功
        if(file.exists(fastqFile)){
          state <- "success"
        }else{
          state <- "error"
        }
      }else{
        # split-3：针对pairend 测序结果，将原文件转换成3个fastq文件，文件名其中一个是*_1.fastq,一个是*_2.fastq,这两个文件里面的reads是成对的，而*.fastq里面的read是剩余不成对的
        t.command <- paste(SRA.dump.fastq, "-s", t.files[i], "-t", threads, "--split-3 --gzip -O", outDir, sep = " ")
        system(t.command)
        fastqFile1 <- paste0(sampleNames[i], "_1.fastq.gz")
        fastqFile2 <- paste0(sampleNames[i], "_2.fastq.gz")
        fastqFile1 <- file.path(outDir, fastqFile1)
        fastqFile2 <- file.path(outDir, fastqFile2)
        
        # 通过检查转换的fastq文件是否存在，来判断转换是否成功
        if(file.exists(fastqFile1) && file.exists(fastqFile2)){
          state <- "success"
        }else{
          state <- "error"
        }
      }
      logout <- c(logout, state)
    }
    logout <- data.frame(sample = sampleNames, state = logout)
    write.table(logout, file = sra2fastq_log_file, col.names = F, row.names = F, quote = F, sep = "\t")
  }else{
    #单样本处理
    if(is.null(outDir)){
      outDir <- dirname(filePath)
    }
    if(pairEND){
      t.command <- paste(SRA.dump.fastq, "-s", filePath, "-t", threads, "--split-files --gzip -O", outDir, sep = " ")
    }else{
      t.command <- paste(SRA.dump.fastq, "-s", filePath, "-t", threads, "--gzip -O", outDir, sep = " ")
    }
    system(t.command)
  }
}