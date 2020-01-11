##  TO DO: 批量运行暴力破解流程的外部函数
## 
##  Author: jeason zhao
###################################################################################

#' @include utilities.R filter.R
#' @import purrr
#' @import survival
#' @import dplyr
#' @import future

# source("~/utilities.R")
# source("~/FilterSignature.R")

#' @name ViolentCrack
#' @description 对signature在各套数据集中进行过滤，要求signature在每套数据集中必须logrank test、单多因素cox p中同时显著。切默认必须进行三种过滤
#' @param clinical.list 一个list，每一个元素是一个dataframe，代表一套数据的临床信息（注意：临床信息必须包含Patient_ID, time, event），
#'   并且将第一个元素作为训练集构建risakscore，其余元素作为验证集
#' @param exp.data.list 一个list，每个元素是一个matrix，代表一套数据的表达谱（注意：行名必须是基因，列名必须是样本）并且将第一个元素作
#'   为训练集，其他作为验证集。特别注意：该list中的信息必须与clinical.data中的信息对应，并且样本数目必须数目和ID必须完全一致
#' @param signatures 矩阵，每个元素必须是字符。包含signature组合信息，每行代表一个signature组合，包含组合中的所有基因。（注意：如果提
#'   供行名，推荐名字为每个signature的gene以下划线“_”形式连接起来的字符串【如:TP53_CDKN2A】，若为其他形式，则必须与cutoff的名字对应一致）
#' @param multi.variable 一个list。每个元素是一个字符型向量，包含参与多因素cox分析的其他临床变量的列名信息
#' @param cutoff 是否使用相同的cutoff区分高低风险组。“same”表示与训练集相同；“diff”表示使用各自的中值作为cutoff值
#' @param candidate 数值，表示在多少套数据集中signature必须是有效，筛选出候选的signature集合。默认为NULL，只保留在所有数据集中有效的signature；
#'   当为任意数值n时，表示计算至少在前n套验证集中logrank有效的signature。如：candidate = 2，表示计算至少在2套验证集中logrank有效的signature
#' @param multiCores 逻辑值，是否进行多核运算，默认为TRUE
#' @param autoCores 逻辑值，是否自动分配内核，默认为TRUE,如果为FALSE，将使用所有可用的线程核
#' @param keepMem 数值，预留百分之多少的内存，防止程序内存不够,默认为0.5，范围是0-1
#' @param outDir 文件目录，用于储存产生的txt文件
#' @param shut.down 逻辑值，是否是在程序中断之后重新运行的，默认为FALSE
#' @param shut.num 数值，用于确认程序中断时signature过滤的位置，只有当shut.down参数为TRUE是有效，默认为NULL。如：当函数在过滤101-200的signature时
#'    中断，则可设置shut.num = 101；程序将继续从第101个signature开始运行
#' @param verbose 是否输出提示信息，默认为TRUE
#' @return 在文件目录中以100万为单位返回暴力破解结果
#' @examples
#' 1. 程序正常运行
#' ViolentCrack(clinical.list, exp.data.list, signatures, multi.variable, cutoff = cutoff,
#'              candidate, multiCores = TRUE, autoCores = TRUE, keepMem = 0.5, outDir, verbose = TRUE))
#' 2. 程序在过滤5000001-6000000个signature时中断后重新运行程序
#' system.time(signature3.results <- ViolentCrack(clinical.list, exp.data.list, signatures, multi.variable, cutoff = cutoff,
#'             candidate, multiCores = TRUE, autoCores = TRUE, keepMem = 0.5, outDir, verbose = TRUE, shut.down = T, shut.num = 5000001)) 
ViolentCrack <- function(clinical.list, exp.data.list, signatures, multi.variable, cutoff = c("same", "diff"),
                         candidate = NULL, multiCores = TRUE, autoCores = TRUE, keepMem = 0.5, outDir,
                         shut.down = FALSE, shut.num = NULL, verbose = TRUE) {
  
  if (verbose) {
    cat("Starting the signature filtering process...", sep = "\n")
  }
  
  ## 设置内部的调节参数
  library(furrr)
  options(future.globals.maxSize = 100.00 * 1024 ^ 3)
  atom <- 10^6 # 以10万为单位分割signature矩阵进行循环
  
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## STEP1:检查数据格式，确定关键变量输入是否正确
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  if (verbose) {
    cat("Checking the data...", sep = "\n")
  }
  
  ## 检查clinical.list,exp.data.list,signatures的输入格式
  stopifnot(.is_list(clinical.list), .is_list(exp.data.list))
  stopifnot(is.matrix(signatures), is.character(signatures))
  
  ## candidate必须不小于2
  if (candidate < 2) {
    stop("`candidate` must >= 2")
  }
  
  ## 检查临床数据和表达谱之间的对应关系
  if (length(clinical.list) != length(exp.data.list)) {
    stop("`clinical.list` and `exp.data.list` are not the same length")
  }
  
  ## 要求数据至少包含一套训练集和一套验证集
  if (length(clinical.list) == 1) {
    stop("The length of `clinical.list` and `exp.data.list` must >= 2")
  }
  
  ## 检查multi.variable参数是否符合条件：(1)是一个list;(2)每个元素包含在每套数据集中多因素cox的协变量的列名;(3)所有变量必须包含在临床数据中
  if (.is_list(multi.variable)) {
    if (all(unlist(purrr::map(multi.variable, is.character)))) {
      if (length(multi.variable) == length(clinical.list)) {
        tmp <- purrr::map2(.x = multi.variable, .y = clinical.list, ~ is.element(.x, colnames(.y)))
        if (!all(unlist(tmp))) {
          stop("The elements of `multi.variable` don't exactly match with `clinical.list`")
        }
      } else {
        stop("`multi.variable` and `clinical.list` have different length")
      }
    } else {
      stop("The elements of `multi.variable` are not all characters")
    }
  } else {
    stop("The argument `multi.variable` must be a list")
  }
  
  ## 确认输出文件夹是否存在。若不存在，则停止函数
  if (!dir.exists(outDir)) {
    stop(outDir, " path doesn't exist, please recheck...")
  }
  
  if (verbose) {
    cat("Adding the rownames for signatures ...", sep = "\n")
  }
  ## 给signature矩阵添加行名
  rownames(signatures) <- apply(signatures, 1, paste, collapse = "_")
  
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## STEP2: 处理临床和表达数据，对样本取交集，便于后续处理
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  if (verbose) {
    cat("Reconstructing the clinical data and gene expression profile ...", sep = "\n")
  }
  ## 对临床样本和表达数据样本取交集，构建新的表达谱和临床数据
  co.patient <- purrr::map2(clinical.list, exp.data.list, ~ intersect(.x$`Patient_ID`, colnames(.y)))
  clinical.list <- purrr::map2(clinical.list, co.patient, ~ dplyr::filter(.x, is.element(Patient_ID, .y)))
  exp.data.list <- purrr::map2(exp.data.list, co.patient, ~ .x[, .y])
  
  ## 使用future包添加并行计算，只需一次声明即可
  worker.cores <- workers(multiCores = multiCores, autoCores = autoCores, keepMem = keepMem)
  future::plan(multicore, workers = worker.cores)
  
  ## 确定是否使用并行计算
  if (multiCores) {
    if (!autoCores) {
      warning("The argument `autoCores` is FALSE, we used all available cores")
    }
  } else {
    warning("The argument `multiCores` is FALSE, all processing will be done in current R session")
  }
  
  ## 以100万为单位，分割signature矩阵，每一份循环进行暴力破解
  tmp <- seq(nrow(signatures))
  m <- ceiling(tmp/atom)
  
  ## 如果程序突然停止，可以根据verbose输出的提示找到停止的点继续执行
  if (shut.down) {
    if (is.null(shut.num)) {
      stop("Argument `shutdown` is TRUE, the `shut.num` can't be NULL.")
    } else {
      k <- ceiling(shut.num/atom)
      for(i in k:max(m)) {
        select.row <- which(m == i)
        if (verbose) {
          cat(paste0("Start filtering the ", min(select.row), "-", max(select.row), " signatures..."), sep = "\n")
        }
        new.outDir <- paste0(outDir, "/", min(select.row), "-", max(select.row)) # 每个单位单独产生一个文件夹用于储存结果
        if (!dir.exists(new.outDir)) {
          dir.create(new.outDir)
        }
        sub.signatures <- signatures[select.row,]
        ## 加入tryCatch，当过滤不出来siganture时自动停止，并且进入下一个循环
        tryCatch(FilterSignature(clinical.list = clinical.list, exp.data.list = exp.data.list, signatures = sub.signatures,
                                 multi.variable = multi.variable, cutoff = cutoff, candidate = candidate, multiCores = multiCores,
                                 autoCores = autoCores, keepMem = keepMem, outDir = new.outDir, verbose = verbose),
                 error = function(e) {cat(paste0("Stop when filtering ", min(select.row), "-", max(select.row), " signatures!"), sep = "\n")})
        gc()
      }
    }
  } else {
    for(i in unique(m)) {
      select.row <- which(m == i)
      if (verbose) {
        cat(paste0("Start filtering the ", min(select.row), "-", max(select.row), " signatures..."), sep = "\n")
      }
      new.outDir <- paste0(outDir, "/", min(select.row), "-", max(select.row)) # 每个单位单独产生一个文件夹用于储存结果
      if (!dir.exists(new.outDir)) {
        dir.create(new.outDir)
      }
      sub.signatures <- signatures[select.row,]
      ## 加入tryCatch，当过滤不出来siganture时自动停止，并且进入下一个循环
      tryCatch(FilterSignature(clinical.list = clinical.list, exp.data.list = exp.data.list, signatures = sub.signatures,
                               multi.variable = multi.variable, cutoff = cutoff, candidate = candidate, multiCores = multiCores,
                               autoCores = autoCores, keepMem = keepMem, outDir = new.outDir, verbose = verbose),
               error = function(e) {cat(paste0("Stop when filtering ", min(select.row), "-", max(select.row), " signatures!"), sep = "\n")})
      gc()
    }
  }
}
