##   TO DO: 暴力破解程序的主体函数，负责每个signature在数据集中的过滤
## 
##   Author: jeason zhao
###################################################################################

#' @include utilities.R
#' @import purrr
#' @import survival
#' @import dplyr
#' @import future

# source("~/utilities.R")

#' @name FilterSignature
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
#' @param verbose 是否输出一系列提示信息，默认为TRUE
#' @return 在文件目录中返回signature的txt文件。每行代表一个signature组合，以tab键为分隔符
#' @examples


FilterSignature <- function(clinical.list, exp.data.list, signatures, multi.variable, cutoff = c("same", "diff"),
                            candidate = NULL, multiCores = TRUE, autoCores = TRUE, keepMem = 0.5, outDir, verbose = TRUE) {
  
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## STEP1: 基于单因素cox系数和基因表达值构建得分模型，区分高低风险组
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  ## 计算signature矩阵中包含的所有基因在训练集中对应的单因素beta系数
  if (verbose) {
    cat("Calculating the beta value ...", sep = "\n")
  }
  beta <- BetaCoefficient(clinical = clinical.list[[1]], exp.data = exp.data.list[[1]], signatures = signatures)
  
  ## 计算risk score = beta1*exp(gene1) + beta2*exp(gene2) + ... + betan*exp(genen)
  if (verbose) {
    cat("Calculating the risk score ...", sep = "\n")
  }
  risk.score <- furrr::future_map(.x = exp.data.list, .f = RiskScore, beta = beta, signatures = signatures)
  
  ## 计算每个数据集的cutoff值：(1)使用训练集的中值；(2)使用各自数据集的中值
  cutoff.list <- vector("list", length = length(clinical.list))
  if (cutoff == "same") {
    cutoff.list[[1]] <- apply(risk.score[[1]], MARGIN = 1, median)
    cutoff.list <- purrr::map(cutoff.list, .f = function(x) x = cutoff.list[[1]])
  } else if (cutoff == "diff") {
    cutoff.list <- purrr::map(risk.score, ~ apply(.x, median, MARGIN = 1))
  }
  
  ## 依据cutoff，区分样本为high_risk和low_risk
  if (verbose) {
    cat("Dividing the sample into high_risk and low_risk groups", sep = "\n")
  }
  label.list <- furrr::future_map2(.x = risk.score, .y = cutoff.list, .f = .label.sample.matrix)
  
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## STEP2: 暴力破解流程，过滤顺序发生修改：要求signatrue必须在训练集
  ## 和前两套验证集中logrank、单因素cox显著，然后重新进行多因素cox过滤
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  if (verbose) {
    cat("Start profiltering ...", sep = "\n")
  }
  
  ##################################################
  ## 前三套数据集中过滤signature，logrank、单因素cox
  #  sub.clinical <- clinical.list[1:3]
  #  sub.exp.data <- exp.data.list[1:3]
  
  ## 提取signature名字
  index <- rownames(label.list[[1]])
  
  ## 循环遍历过滤
  for (j in 1:length(clinical.list[1:3])) {
    ## 提取相应临床信息和label信息
    clinical <- clinical.list[[j]]
    label.sample <- label.list[[j]][index,]
    
    #################################################
    ## logrank 过滤
    #################################################
    
    ## 如果取出的signature对应的label只有一行，将其转化成矩阵形式
    if (is.vector(label.sample)) {
      label.sample <- label.sample[match(names(label.sample), clinical[, "Patient_ID"])]
      label.sample <- t(as.matrix(label.sample))
      rownames(label.sample) <- index
    } else if (is.matrix(label.sample)) {
      label.sample <- label.sample[, match(colnames(label.sample), clinical[, "Patient_ID"])]
    }
    
    ## 对所有signature进行logrank test，不能计算p值的将p值设为NA
    if (verbose) {
      if (j == 1) {
        cat("Start the logrank test in training dataset ...", sep = "\n")
      } else {
        temp <- sprintf("Starting the logrank test in the %i validation data set ...", j-1)
        cat(temp, sep = "\n")
      }
    }
    logrank.p <- do.call(c, furrr::future_map(.x = split(label.sample, row(label.sample)), function(x) {
      fit <- try(p <- logRank(data = clinical, time = clinical$time,
                              event = all.clinical$event, variate = x), silent = TRUE)
      if ('try-error' %in% class(fit)) {
        p <- NA
      } else {
        p <- p
      }
      return(p)
    }))
    
    ## 判断每个signature的logrank test的p值是否显著，过滤掉不显著的signature
    index <- rownames(label.sample)[.compare_na(logrank.p < 0.05)]
    if (verbose) {
      cat(paste0("left ", length(index), " signatures"), sep = "\n")
    }
    
    ## 判断在过滤p值不显著的signature后是否仍有剩余。若没有剩余的signature，则跳出程序
    if (is.character0(index)) {
      if (j == 1) {
        cat("When the logrank test is performed in the training data set, no signature meets the condition.", sep = "\n")
      } else {
        tmp <- sprintf("When the logrank test is performed in the %i validation data set, no signature meets the condition.", j-1)
        cat(tmp, sep = "\n")
      }
      break # 跳出循环
    } else {
      ## 过滤label.list，将不符合条件的signature相关信息过滤掉
      label.list <- purrr::map(label.list, ~ .x[index,])
      if (length(index) < 2) {
        label.list <- purrr::map(label.list, ~ t(as.matrix(.x)))
        label.list <- purrr::map(label.list, .f = function(x) {
          rownames(x) <- index
          x
        })
      }
    }
    
    #################################################
    ## 单因素cox过滤
    #################################################
    
    ## 依据过滤后的signature，提取label信息。如果取出的signature对应的label只有一行，将其转化成矩阵形式
    label.sample <- label.list[[j]][index,]
    if (is.vector(label.sample)) {  
      label.sample <- label.sample[match(names(label.sample), clinical[, "Patient_ID"])]
      label.sample <- t(as.matrix(label.sample))
      rownames(label.sample) <- index
    } else if(is.matrix(label.sample)) {
      label.sample <- label.sample[, match(colnames(label.sample), clinical[, "Patient_ID"])]
    }
    
    ## 对所有signature进行单因素cox分析,不能计算p值的将p值设为NA
    if (verbose) {
      if (j == 1) {
        cat("Start the univariate Cox analysis in training dataset ...", sep = "\n")
      } else {
        temp <- sprintf("Starting the univariate Cox analysis in the %i validation data set ...", j-1)
        cat(temp, sep = "\n")
      }
    }
    uni.p <- do.call(c, furrr::future_map(.x = split(label.sample, row(label.sample)), function(x) {
      fit <- try(p <- uniVariatePvalue(data = clinical, time = clinical$time,
                                       event = clinical$event, variate = x)["Pr(>|z|)"], silent = TRUE)
      if ('try-error' %in% class(fit)) {
        p <- NA
      } else {
        p <- p
      }
      return(p)
    }))
    
    ## 判断每个signature的单因素cox的p值是否显著，过滤掉不显著的signature
    index <- rownames(label.sample)[.compare_na(uni.p < 0.05)]
    if (verbose) {
      cat(paste0("left ", length(index), " signatures"), sep = "\n")
    }
    ## 判断在过滤p值不显著的signature后是否仍有剩余。若没有剩余的signature，则跳出程序
    if (is.character0(index)) {
      if (j == 1) {
        cat("When the univariate Cox analysis is performed in the training data set, no signature meets the condition.", sep ="\n")
      } else {
        tmp <- sprintf("When the univariate Cox analysis is performed in the %i validation data set, no signature meets the condition.", j-1)
        cat(tmp, sep = "\n")
      }
      break # 跳出循环
    } else {
      ## 过滤label.list，将不符合条件的signature相关信息过滤掉
      label.list <- purrr::map(label.list, ~ .x[index,])
      if (length(index) < 2) {
        label.list <- purrr::map(label.list, ~ t(as.matrix(.x)))
        label.list <- purrr::map(label.list, .f = function(x) {
          rownames(x) <- index
          x
        })
      }
    }
  }
  
  if (!is.character0(index)) {
    index <- rownames(label.list[[1]])
    ## 开始重新进行多因素cox计算（在训练集和前两套验证集中）
    for (j in 1:length(clinical.list[1:3])) {
      ## 提取相应临床信息和label信息
      clinical <- clinical.list[[j]]
      label.sample <- label.list[[j]][index,]
      
      if (is.vector(label.sample)) {  
        label.sample <- label.sample[match(names(label.sample), clinical[, "Patient_ID"])]
        label.sample <- t(as.matrix(label.sample))
        rownames(label.sample) <- index
      } else if(is.matrix(label.sample)) {
        label.sample <- label.sample[, match(colnames(label.sample), clinical[, "Patient_ID"])]
      }
      
      ## 提取多因素cox需要的其他临床协变量名称 
      if (verbose) {
        if (j == 1) {
          cat("Start the multivariate Cox analysis in training dataset ...", sep = "\n")
        } else {
          temp <- sprintf("Starting the multivariate Cox analysis in the %i validation data set ...", j-1)
          cat(temp, sep = "\n")
        }
      }
      multi.var <- multi.variable[[j]]
      
      ## 对所有signature进行多因素cox分析,不能计算p值的将p值设为NA
      tmp <- tibble::rownames_to_column(as.data.frame(t(label.sample)), var = "Patient_ID")
      all.clinical <- dplyr::left_join(clinical, tmp, by = "Patient_ID")
      multi.p <- do.call(c, furrr::future_map(rownames(label.sample), function(x) {
        fit <- try(p <- multiVariatePvalue(data = all.clinical, time = all.clinical$time, 
                                           event = all.clinical$event, 
                                           variate = c(x, multi.var))["Pr(>|z|)"], silent = TRUE)
        if ('try-error' %in% class(fit)) {
          p <- NA
        } else {
          p <- p
        }
        return(p)
      }))
      
      ## 判断每个signature的多因素cox的p值是否显著，过滤掉不显著的signature
      index <- rownames(label.sample)[.compare_na(multi.p < 0.05)]
      if (verbose) {
        cat(paste0("left ", length(index), " signatures"), sep = "\n")
      }
      
      ## 判断在过滤p值不显著的signature后是否仍有剩余。若没有剩余的signature，则跳出程序
      if (is.character0(index)) {
        if (j == 1) {
          cat("When the multivariate Cox analysis is performed in the training data set, no signature meets the condition.", sep = "\n")
        } else {
          tmp <- sprintf("When the multivariate Cox analysis is performed in the %i validation data set, no signature meets the condition.", j-1)
          cat(tmp, sep = "\n")
        }
        break # 跳出循环
      } else {
        ## 根据candidate参数，确定是否计算输出候选signature
        if (!is.null(candidate)) {
          if (!(j < candidate + 1)) {
            candidate.sig <- do.call(rbind, strsplit(index, split = "_"))
            write.table(candidate.sig, file = paste0(outDir, "/candidate", j - candidate, ".txt"),
                        sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
          }
        }
        
        ## 过滤label.list，将不符合条件的signature相关信息过滤掉
        label.list <- purrr::map(label.list, ~ .x[index,])
        if (length(index) < 2) {
          label.list <- purrr::map(label.list, ~ t(as.matrix(.x)))
          label.list <- purrr::map(label.list, .f = function(x) {
            rownames(x) <- index
            x
          })
        }
      }
    }
    
    if (!is.character0(index)) {
      if (verbose) {
        cat("profiltering finished ...", sep = "\n")
        cat("Start the filtering ...", sep = "\n")
      }
      index <- rownames(label.list[[1]])
      ## 开始后续的正常过滤流程
      for (j in 4:length(clinical.list)) {
        if (verbose) {
          cat(paste0("Satrting the process in ", j - 1, " validation data set ..."), sep = "\n")
        }
        ## 提取相应临床信息和label信息
        clinical <- clinical.list[[j]]
        label.sample <- label.list[[j]][index,]
        
        #################################################
        ## logrank 过滤
        #################################################
        
        ## 如果取出的signature对应的label只有一行，将其转化成矩阵形式
        if (is.vector(label.sample)) {
          label.sample <- label.sample[match(names(label.sample), clinical[, "Patient_ID"])]
          label.sample <- t(as.matrix(label.sample))
          rownames(label.sample) <- index
        } else if (is.matrix(label.sample)) {
          label.sample <- label.sample[, match(colnames(label.sample), clinical[, "Patient_ID"])]
        }
        
        ## 对所有signature进行logrank test，不能计算p值的将p值设为NA
        if (verbose) {
          cat("Starting the logrank test ...", sep = "\n")
        }
        logrank.p <- do.call(c, furrr::future_map(.x = split(label.sample, row(label.sample)), function(x) {
          fit <- try(p <- logRank(data = clinical, time = clinical$time,
                                  event = all.clinical$event, variate = x), silent = TRUE)
          if ('try-error' %in% class(fit)) {
            p <- NA
          } else {
            p <- p
          }
          return(p)
        }))
        
        ## 判断每个signature的logrank test的p值是否显著，过滤掉不显著的signature
        index <- rownames(label.sample)[.compare_na(logrank.p < 0.05)]
        if (verbose) {
          cat(paste0("left ", length(index), " signatures"), sep = "\n")
        }
        
        ## 判断在过滤p值不显著的signature后是否仍有剩余。若没有剩余的signature，则跳出程序
        if (is.character0(index)) {
          tmp <- sprintf("When the logrank test is performed in the %i validation data set, no signature meets the condition.", j-1)
          cat(tmp, sep = "\n")
          break # 跳出循环
        } else {
          ## 过滤label.list，将不符合条件的signature相关信息过滤掉
          label.list <- purrr::map(label.list, ~ .x[index,])
          if (length(index) < 2) {
            label.list <- purrr::map(label.list, ~ t(as.matrix(.x)))
            label.list <- purrr::map(label.list, .f = function(x) {
              rownames(x) <- index
              x
            })
          }
          
          #################################################
          ## 单因素cox过滤
          #################################################
          
          ## 依据过滤后的signature，提取label信息。如果取出的signature对应的label只有一行，将其转化成矩阵形式
          label.sample <- label.list[[j]][index,]
          if (is.vector(label.sample)) {  
            label.sample <- label.sample[match(names(label.sample), clinical[, "Patient_ID"])]
            label.sample <- t(as.matrix(label.sample))
            rownames(label.sample) <- index
          } else if(is.matrix(label.sample)) {
            label.sample <- label.sample[, match(colnames(label.sample), clinical[, "Patient_ID"])]
          }
          
          ## 对所有signature进行单因素cox分析,不能计算p值的将p值设为NA
          if (verbose) {
            cat("Starting the univariate Cox analysis ...", sep = "\n")
          }
          uni.p <- do.call(c, furrr::future_map(.x = split(label.sample, row(label.sample)), function(x) {
            fit <- try(p <- uniVariatePvalue(data = clinical, time = clinical$time,
                                             event = clinical$event, variate = x)["Pr(>|z|)"], silent = TRUE)
            if ('try-error' %in% class(fit)) {
              p <- NA
            } else {
              p <- p
            }
            return(p)
          }))
          
          ## 判断每个signature的单因素cox的p值是否显著，过滤掉不显著的signature
          index <- rownames(label.sample)[.compare_na(uni.p < 0.05)]
          if (verbose) {
            cat(paste0("left ", length(index), " signatures"), sep = "\n")
          }
          
          ## 判断在过滤p值不显著的signature后是否仍有剩余。若没有剩余的signature，则跳出程序
          if (is.character0(index)) {
            tmp <- sprintf("When the univariate Cox analysis is performed in the %i validation data set, no signature meets the condition.", j-1)
            cat(tmp, sep = "\n")
            break # 跳出循环
          } else {
            ## 过滤label.list，将不符合条件的signature相关信息过滤掉
            label.list <- purrr::map(label.list, ~ .x[index,])
            if (length(index) < 2) {
              label.list <- purrr::map(label.list, ~ t(as.matrix(.x)))
              label.list <- purrr::map(label.list, .f = function(x) {
                rownames(x) <- index
                x
              })
            }
            
            #################################################
            ## 多因素cox过滤
            #################################################
            
            ## 依据过滤后的signature，提取label信息。如果取出的signature对应的label只有一行，将其转化成矩阵形式
            label.sample <- label.list[[j]][index,]
            if (is.vector(label.sample)) {  
              label.sample <- label.sample[match(names(label.sample), clinical[, "Patient_ID"])]
              label.sample <- t(as.matrix(label.sample))
              rownames(label.sample) <- index
            } else if(is.matrix(label.sample)) {
              label.sample <- label.sample[, match(colnames(label.sample), clinical[, "Patient_ID"])]
            }
            
            ## 提取多因素cox需要的其他临床协变量名称 
            if (verbose) {
              cat("Starting the multivariate Cox analysis ...", sep = "\n")
            }
            multi.var <- multi.variable[[j]]
            
            ## 对所有signature进行多因素cox分析,不能计算p值的将p值设为NA
            tmp <- tibble::rownames_to_column(as.data.frame(t(label.sample)), var = "Patient_ID")
            all.clinical <- dplyr::left_join(clinical, tmp, by = "Patient_ID")
            multi.p <- do.call(c, furrr::future_map(rownames(label.sample), function(x) {
              fit <- try(p <- multiVariatePvalue(data = all.clinical, time = all.clinical$time, 
                                                 event = all.clinical$event, 
                                                 variate = c(x, multi.var))["Pr(>|z|)"], silent = TRUE)
              if ('try-error' %in% class(fit)) {
                p <- NA
              } else {
                p <- p
              }
              return(p)
            }))
            
            ## 判断每个signature的多因素cox的p值是否显著，过滤掉不显著的signature
            index <- rownames(label.sample)[.compare_na(multi.p < 0.05)]
            if (verbose) {
              cat(paste0("left ", length(index), " signatures"), sep = "\n")
            }
            
            ## 根据candidate参数，确定是否计算输出候选signature
            if (is.null(candidate)) {
              if (is.character0(index)) {
                cat("There was no signature is filtered out!", sep = "\n")
                break # 跳出循环
              }
            } else {
              if (is.character0(index)) {
                if (j > candidate +1) {
                  tmp <- sprintf("No signature meets condition after the %i validation data set.", j-1)
                  cat(tmp, sep = "\n")
                  break # 跳出循环
                } else {
                  cat("There was no signature is filtered out!", sep = "\n")
                  break # 跳出循环
                }
              } else {
                if (!(j < candidate + 1)) {
                  candidate.sig <- do.call(rbind, strsplit(index, split = "_"))
                  write.table(candidate.sig, file = paste0(outDir, "/candidate", j - candidate, ".txt"),
                              sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
                }
              }
            }
          }
        }
      }
    }
  }
}

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##help function
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Check if is a list
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
.is_list <- function(x){
  inherits(x, "list")
}

# Check if a character is 0
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
is.character0 <- function(x) {
  is.character(x) && length(x) == 0L
}

# convert NA to FALSE when compare
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
.compare_na <- function(expr, type = FALSE) {
  expr[is.na(expr)] <- type
  expr
}

# 区分高低风险的矩阵操作
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
.label.sample.matrix <- function(risk.score.mat, cutoff) {
  risk.score.list <- split(risk.score.mat, row(risk.score.mat))
  label.sample <- do.call(rbind, purrr::map2(.x = risk.score.list, 
                                             .y = cutoff, .f = LableSample))
  dimnames(label.sample) <- dimnames(risk.score.mat)
  label.sample
}