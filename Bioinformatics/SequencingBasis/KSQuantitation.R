####################################################################################
# 整体功能简介：
# 使用kallisto/salmon进行转录组定量表达
###################################################################################
# 重要使用函数简介
# 代码中共包括两个函数：
# 1）make.index，用于建立kallisto/salmon的index文件
# 2）ks.quantitation，实现kallisto/salmon的定量转录本表达

###################################################################################
#函数一：make.index
#' @@ 建立kallisto或salmon要求的index文件
#' @param transcriptomeFile 参考转录组文件绝对地址，fq或fq.gz等gz压缩格式均可，暂不支持其他压缩形式。推荐使用genecode，打开 https://www.gencodegenes.org/human/ 网址，选择Fasta files--> Transcript sequences下载
#' salmon和kallisto采用的是参考转录组，不建议使用参考基因组。如"/pub6/temp/xlw/kallisto_linux-v0.44.0/gencode.v29.transcripts.fa.gz"。
#' @param outDir index文件输出目录，需要事先建好。
#' @param type 软件使用类型，当前支持"kallisto"和"salmon"，注意输入格式要保持完全正确。
#' @param k 有效匹配的最小可接受长度，当前默认为31。注意这个参数需要用户确认read的长度范围，当fastq文件中大量read长度小于31时，软件处理会出现问题。
#' kallisto软件默认为31，也是支持的最大值。
#' salmon官方文档说明，k设置为31对于read长度为75bp或者更长效果要好，如果用户需要处理更短的read，可以考虑设置一个更小的k值，但并没有给出具体设置方案。

#' @return index文件。kallisto index是一个idx文件结尾的文件，而salmon index是一个文件夹

#' @note 
#' 当前使用Salmon-v0.9.1，kallisto 0.44.0版本
#' 运行前需要确认这两个软件是否存在在给定的运行路径下:salmon.path = "/pub5/xiaoyun/BioSoftware/Salmon-v0.9.1/bin/salmon";kallisto.path = "/pub5/xiaoyun/BioSoftware/kallisto_linux-v0.44.0/kallisto"
#' 对于salmon，支持两种方式建立index（Quasi-index and FMD-index-based modes），当前使用更优的Quasi-index模式
#' salmon建立index，可能会出现大片黄色警告，一般情况都可以忽略：通常是表示有些转录本长度不满足31的要求；salmon还会删除定义为重复的转录本，记录在index文件夹下的duplicate_clusters.tsv

#' example
#' make.index(transcriptomeFile = "/pub6/temp/xlw/kallisto_linux-v0.44.0/gencode.v29.transcripts.fa.gz", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test", type = "salmon", k = 31)
#' make.index(transcriptomeFile = "/pub6/temp/xlw/kallisto_linux-v0.44.0/gencode.v29.transcripts.fa.gz", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test", type = "kallisto", k = 31)

make.index <- function(transcriptomeFile, outDir, type = c("kallisto", "salmon"), k = 31){
  
  salmon.path <- "/pub5/xiaoyun/BioSoftware/Salmon-v0.9.1/bin/salmon"
  kallisto.path <- "/pub5/xiaoyun/BioSoftware/kallisto_linux-v0.44.0/kallisto"
  
  if(!dir.exists(outDir)){
    stop(outDir, "目录不存在，请确认...")
  }
  
  # 获取参考转录本名字用于命名index文件
  ##格式通常为fa/fa.gz/fasta/fasta.gz
  gz <- grep("gz$", transcriptomeFile)
  if(length(gz)>0){
    faName <- gsub("\\.f[^.]*\\.gz$", "", transcriptomeFile, ignore.case = T)
  }else{
    faName <- gsub("\\.f[^.]*$", "", transcriptomeFile, ignore.case = T)
  }
  faName <- gsub(".*/", "", faName)
  
  if(type == "kallisto"){
    outIndexFile <- paste0(faName, ".idx")
    outIndexFile <- file.path(outDir, outIndexFile)
    index.command <- paste(kallisto.path, "index", "-i", outIndexFile, "-k", k, transcriptomeFile, sep = " ")
  }
  if(type == "salmon"){
    outIndexFile <- faName
    outIndexFile <- file.path(outDir, outIndexFile)
    index.command <- paste(salmon.path, "index -t", transcriptomeFile, "-i", outIndexFile, "--type quasi", "-k", k, sep = " ")
  }
  
  print(paste0(type, " index building start..."))
  system(index.command)   
}

###################################################################################
#函数二：ks.quantitation
#' @param fastqDir fastq文件目录，如/IData/CancerOMICS/SkinCancer/Riaz_Cell_2017/RNA_Seq/trim_galore
#' @param is.TrimGalore TRUE/FALSE，是否是经过trim galore处理过的fastq数据
#' @param pattern fastq文件的结尾模式，用于匹配fastqDir目录下的fastq文件，如"fastq/fq.gz"或"fastq/fq"，不能为".fastq/fq.gz"或".fastq/fq"
#' @param type 软件使用类型，当前支持kallisto和salmon。注意输入要保持正确
#' @param index.path 转录组index文件或目录，
#' salmon要求index目录，如/home/longzl2017/dataScience/PipeLine/RNA-seq/test/gencode_v28_transcripts_index
#' kallisto要求index文件，如/home/longzl2017/dataScience/PipeLine/RNA-seq/test/Homo_gencodeV29_transcripts.idx
#' @param outDir 输出目录，要求具有可写权限，如/IData/CancerOMICS/SkinCancer/Riaz_Cell_2017/RNA_Seq/quantitation_result。如果不存在，则会自动建立
#' @param pairEND: TRUE/FALSE, 是否是双末端测序, 这个一定要仔细查证之后再给定，因为kallisto或salmon处理单双末端参数不同
#' @param fldMean 估计的片段平均长度(mean fragment lenth of the sequencing library)，只适用单末端模式，且是kallisto必须指定的重要参数，会直接影响定量的准确性。
#' 当前采用kallisto官网的默认值100，实际使用时最好进一步查证，根据自己数据的实际情况进行修改设定。软件fastqc可以查看read长度分布情况
#' salmon官网上并没有强制要求单末端时必须指定该参数，所以type参数为salmon时，允许设置为NULL，但是明确说明了该参数很重要，会直接影响定量的准确性，使不使需要自己斟酌
#' @param fldSD 估计的片段长度的标准差(standard deviation of the fragment lenth distribution of the sequencing library)，只适用单末端模式，且是kallisto必须指定的重要参数，会直接影响定量的准确性。
#' salmon官网上并没有强制要求单末端时必须指定该参数，所以type参数为salmon时，允许设置为NULL，但是明确说明了该参数很重要，会直接影响定量的准确性，使不使需要自己斟酌
#' @param strand.specific 是否使用链特异模式运行kallisto，可选参数，默认为 "" ;支持"fr-stranded"或"rf-stranded"输入，"fr-stranded"指用链特异模式运行kallisto
#' "rf-stranded" 也指链特异，但first read 比对到转录本的反义链上。只支持kallisto
#' @param is.Bias 是否进行偏性矫正，默认为FALSE，不进行偏性矫正。kallisto支持基于sequence的bias矫正；salmon支持矫正序列特异的bias和fragment-level GC的bias。个人推荐使用，用户可以进一步查证自行选择。
#' @param is.BAMoutput 是否需要BAM格式的输出，默认FALSE，不输出BAM文件，只支持kallisto
#' @param extraParameter kallisto或salmon的自定义扩展参数字符串，可以加入不在函数设定参数内的运行参数，默认为NULL。输入需要完全按照参数组合方式排列，如"--numBootstraps 40"
#' @param threads 线程数目，默认为20 

#' @return 转录本定量结果文件。
#### kallisto会为每个样本3个文件： 详细说明见 --- http://210.46.85.145/showtopic-2243.aspx
#' abundances.tsv 是转录本表达丰度定量文件，不包含bootstrap estimates. 使用--plaintext 参数输出文本格式的abundance estimates。或者用h5dump命令把输出的HDF5文件转成文本格式。第一行是每列的列名，包含estimated counts, TPM, effective length.
#' run_info.json 是一个 json 文件，包含程序运行参数log的信息。

#### salmon会为每个样本产生一个quant.sf(定量结果)和很多其他记录文件：详细说明见 --- http://210.46.85.145/showtopic-2272.aspx
#' quant.sf 是转录本表达丰度定量文件，包含Name(转录本名字)、Length(length of the target transcript in nucleotides)、EffectiveLength(估计转录本有效长度)、TPM(估计的转录本TPM表达丰度)、NumReads(估计的比对到转录本上的read数目)
#' cmd_info.json 是程序运行参数记录文件
#' lib_format_counts.json记录了salmon估计的样本中所有read的文库协议分布情况
#' 其他文件可从 https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats 进一步查阅

#' @note 
#' 内嵌的软件执行路径： 运行前需要确认一下
#' Salmon: /pub5/xiaoyun/BioSoftware/Salmon-v0.9.1/bin/salmon
#' kallisto: /pub5/xiaoyun/BioSoftware/kallisto_linux-v0.44.0/kallisto
#' Samtools: /pub5/xiaoyun/BioSoftware/Samtools-1.3/bin
#' 
#' 当前使用Salmon-v0.9.1，kallisto 0.44.0版本。salmon相对比kallisto运行速度更快
#' kallisto和salmon定量结果是转录本水平的，所以我们需要额外的一步转换(transcript2gene.quant.R)，将转录本表达转换为基因水平的表达
#' 对于多个run的样本，需要在定量前将同一个样本多个run合并
#' 
#### salmon
#' salmon定量默认使用参数：-l A 自动确定文库类型；
#' 使用salmon定量时输出控制台会显示出一段黄色输出：显示了有多少比例的read小于设定的K值和多少read比对上了，这段信息对于参数是否设置正确具有很好的指导意义，需要着重注意!!!
#### kallisto
#' kallisto定量默认使用参数：--plaintext 使用普通文本替代HDF5文件输出

#' example
#' single-end kallisto
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq", is.TrimGalore = F, pattern = "fastq.gz", type = "kallisto", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/Homo_gencodeV29_transcripts.idx", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq/kresult", pairEND = F, fldMean = 100, fldSD = 20)
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/trimGalore_result", is.TrimGalore = T, type = "kallisto", pattern = "fq.gz$", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/Homo_gencodeV29_transcripts.idx", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/trimGalore_result/kresult", pairEND = F, fldMean = 100, fldSD = 20)
#' 
#' single-end salmon
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq", is.TrimGalore = F, pattern = "fastq.gz", type = "salmon", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/gencode_v29_transcripts_index", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/fastq/sresult", pairEND = F, fldMean = 100, fldSD = 20)
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/trimGalore_result", is.TrimGalore = T, type = "salmon", pattern = "fq.gz$", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/gencode_v29_transcripts_index", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/single/trimGalore_result/sresult", pairEND = F, fldMean = 100, fldSD = 20)
#' 
#' paired-end kallisto
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/fastq", is.TrimGalore = F, type = "kallisto", pattern = "fastq.gz", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/Homo_gencodeV29_transcripts.idx", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/fastq/kresult", pairEND = T)
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/trimGalore_result", is.TrimGalore = T, type = "kallisto", pattern = "fq.gz$", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/Homo_gencodeV29_transcripts.idx", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/trimGalore_result/kresult", pairEND = T)
#' 
#' paired-end salmon
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/fastq", is.TrimGalore = F, type = "salmon", pattern = "fastq.gz", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/gencode_v29_transcripts_index", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/fastq/sresult", pairEND = T)
#' ks.quantitation(fastqDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/trimGalore_result", is.TrimGalore = T, type = "salmon", pattern = "fq.gz$", index.path = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/gencode_v29_transcripts_index", outDir = "/home/longzl2017/dataScience/PipeLine/RNA-seq/test/paired/trimGalore_result/sresult", pairEND = T)

ks.quantitation <- function(fastqDir, is.TrimGalore = FALSE, pattern = "fq.gz$", type = c("kallisto", "salmon"), index.path, outDir, pairEND, fldMean = 100, fldSD = 20, is.Bias = FALSE, is.BAMoutput = FALSE, strand.specific = "", extraParameter = NULL, threads = 20){
  
  # 软件运行路径
  salmon.path <- "/pub5/xiaoyun/BioSoftware/Salmon-v0.9.1/bin/salmon"
  kallisto.path <- "/pub5/xiaoyun/BioSoftware/kallisto_linux-v0.44.0/kallisto"
  samtools.path <- "/pub5/xiaoyun/BioSoftware/Samtools-1.3/bin"
  
  # 检查文件是否存在
  if(!dir.exists(fastqDir)){
    stop(fastqDir, " not exist!please check...")
  }
  
  # 如果outDir不存在且不为空，则自动建立存储目录
  if(!dir.exists(outDir)){
    if(!is.null(outDir)){
      dir.create(outDir)
    }
  }
  
  # 用于指示type参数输出是否有误
  if(!(type %in% c("kallisto", "salmon"))){
    stop("请确认输入的type参数是否正确!")
  }
  
  # 双末端
  if(pairEND){
    paired.fastqfile1 <- list.files(fastqDir, pattern=paste0("_1.", pattern), full.names=TRUE, recursive=F)
    paired.fastqfile2 <- list.files(fastqDir, pattern=paste0("_2.", pattern), full.names=TRUE, recursive=F)
    paired.sampleName <- list.files(fastqDir, pattern=paste0("_1.", pattern), full.names=F, recursive=F)  
    
    # 衔接trim galore处理结果，获得更准确的样本名
    if(is.TrimGalore){
      paired.sampleName <- gsub(paste0("_1_val_1.", pattern), "", paired.sampleName)
    }else{
      paired.sampleName <- gsub(paste0("_1.", pattern), "", paired.sampleName)
    }
    
    if(type == "kallisto"){
      # 检查index是否存在且格式是否正确
      if(!file.exists(index.path) || (length(grep(".idx$", index.path)) == 0)){
        stop(index.path, " 不存在或格式不符合!please check...")
      }
      
      kallisto.command <- paste(kallisto.path, "quant -i", index.path, "-t", threads, "--plaintext", sep = " ")
      
      # 是否进行偏性矫正
      if(is.Bias){
        kallisto.command <- paste0(kallisto.command, " --bias") 
      }
      
      # 是否是链特异
      if(nchar(strand.specific) > 0){
        kallisto.command <- paste0(kallisto.command, " --", strand.specific)
      }
      
      # 加入扩展参数
      if(!is.null(extraParameter)){
        kallisto.command <- paste(kallisto.command, extraParameter, sep = " ") 
      }
      
      for(i in 1:length(paired.sampleName)){
        # 构建输出目录
        sampleDir <- file.path(outDir, paired.sampleName[i])
        if(!dir.exists(sampleDir)){
          dir.create(sampleDir)                                
        }
        
        ##---------如果要输出BAM文件
        if(is.BAMoutput){
          paired.command <- paste0(kallisto.command, " -o ", sampleDir, " --pseudobam ", paired.fastqfile1[i], " ", paired.fastqfile2[i], " | samtools view -Sb - > ", sampleDir, "/out.bam")
        }else{
          paired.command <- paste(kallisto.command, "-o", sampleDir, paired.fastqfile1[i], paired.fastqfile2[i], sep = " ")
        } 
        
        print(paired.command)
        system(paired.command)
      }
    }
    
    if(type == "salmon"){
      # 检查index是否存在且格式是否正确
      if(!file.exists(index.path) || (length(grep(".idx$", index.path)) > 0)){
        stop(index.path, " 不存在或格式不符合!please check...")
      }
      
      salmon.command <- paste(salmon.path, "quant -i", index.path, "-p", threads, "-l A", sep = " ")
      
      # 是否进行偏性矫正
      if(is.Bias){
        salmon.command <- paste0(salmon.command, " --seqBias --gcBias") 
      }
      
      # 加入扩展参数
      if(!is.null(extraParameter)){
        salmon.command <- paste(salmon.command, extraParameter, sep = " ")
      }
      
      for(i in 1:length(paired.sampleName)){
        # 构建输出目录
        sampleDir <- file.path(outDir, paired.sampleName[i])
        if(!dir.exists(sampleDir)){
          dir.create(sampleDir)                                
        }
        paired.command <- paste(salmon.command, "-1", paired.fastqfile1[i], "-2", paired.fastqfile2[i], "-o", sampleDir, sep = " ")
        
        print(paired.command)
        system(paired.command)
      }
    }
    
  }else{
    # 单末端
    single.fastqfile <- list.files(fastqDir, pattern=pattern, full.names=TRUE, recursive=F)
    single.sampleName <- list.files(fastqDir, pattern=pattern, full.names=F, recursive=F)
    
    # 衔接trim galore处理结果，获得更准确的样本名
    if(is.TrimGalore){
      single.sampleName <- gsub(paste0("_trimmed.", pattern), "", single.sampleName) 
    }else{
      single.sampleName <- gsub(paste0(".", pattern), "", single.sampleName) 
    }
    
    # kallisto定量
    if(type == "kallisto"){
      # 检查index是否存在且格式是否正确
      if(!file.exists(index.path) || (length(grep(".idx$", index.path)) == 0)){
        stop(index.path, " 不存在或格式不符合!please check...")
      }
      
      kallisto.command <- paste(kallisto.path, "quant -i", index.path, "-t", threads, "--single -l", fldMean, "-s", fldSD, "--plaintext", sep = " ")
      
      # 是否进行偏性矫正
      if(is.Bias){
        kallisto.command <- paste0(kallisto.command, " --bias") 
      }
      
      # 是否是链特异
      if(nchar(strand.specific) > 0){
        kallisto.command <- paste0(kallisto.command, " --", strand.specific)
      }
      
      # 加入扩展参数
      if(!is.null(extraParameter)){
        kallisto.command <- paste(kallisto.command, extraParameter, sep = " ") 
      }
      
      for(i in 1:length(single.sampleName)){
        # 构建输出目录
        sampleDir <- file.path(outDir, single.sampleName[i])
        if(!dir.exists(sampleDir)){
          dir.create(sampleDir)                                
        }
        ##---------如果要输出BAM文件
        if(is.BAMoutput){
          single.command <- paste0(kallisto.command, " -o ", sampleDir, " --pseudobam ", single.fastqfile[i], " | samtools view -Sb - > ", sampleDir, "/out.bam")
        }else{
          single.command <- paste(kallisto.command, "-o", sampleDir, single.fastqfile[i], sep = " ")
        }
        
        print(single.command)                        
        system(single.command)
      }
    }
    
    # salmon定量
    if(type == "salmon"){
      # 检查index是否存在且格式是否正确
      if(!file.exists(index.path) || (length(grep(".idx$", index.path)) > 0)){
        stop(index.path, " 不存在或格式不符合!please check...")
      }
      
      # 平均片段长度
      if(is.null(fldMean)){
        salmon.command <- paste(salmon.path, "quant -i", index.path, "-p", threads, "-l A", sep = " ")
      }else{
        salmon.command <- paste(salmon.path, "quant -i", index.path, "-p", threads, "-l A --fldMean", fldMean, sep = " ")
      }
      
      # 片段长度方差
      if(!is.null(fldSD)){
        salmon.command <- paste(salmon.command, "--fldSD", fldSD, sep = " ")
      }
      
      # 是否进行偏性矫正
      if(is.Bias){
        salmon.command <- paste0(salmon.command, " --seqBias --gcBias") 
      }
      
      # 加入扩展参数
      if(!is.null(extraParameter)){
        salmon.command <- paste(salmon.command, extraParameter, sep = " ")
      }
      
      for(i in 1:length(single.sampleName)){
        # 构建输出目录
        sampleDir <- file.path(outDir, single.sampleName[i])
        if(!dir.exists(sampleDir)){
          dir.create(sampleDir)                                
        }
        single.command <- paste(salmon.command, "-r", single.fastqfile[i], "-o", sampleDir, sep = " ")
        
        print(single.command)                        
        system(single.command)
      }
    }
  }
}


