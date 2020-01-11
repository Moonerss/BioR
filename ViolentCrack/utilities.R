##  TO DO: 暴力破解程序中使用的基础函数，如logrank p值计算
## 
##  Author: jeason zhao
###################################################################################
#' @import parallel
#' @import survival
#' @import future

#' @name sysmem
#' @description 返回系统总内存数,单位Mb
#' @returnType 
#' @return 
sysmem <- function() {
  t1 <- system("cat /proc/meminfo", intern = TRUE);
  t1 <- gsub("[[:alpha:][:punct:][:space:]]", "", t1[1])
  return(as.numeric(t1) / 1024) ## 单位Mb
}

#' @name memR
#' @description 当前R内可用内存为多少,单位Mb
#' @returnType 
#' @return 
memR <- function() {
  bit <- 8L * .Machine$sizeof.pointer
  if (!(bit == 32L || bit == 64L)) {
    stop("Unknown architecture", call. = FALSE)
  }
  
  node_size <- if (bit == 32L) 28L else 56L
  
  usage <- gc()
  sum(usage[, 1] * c(node_size, 8)) / (1024 ^ 2) ## 单位Mb
}

#' @name memLinux
#' @description 返回当前R在linux中占用的内存，单位Mb
#' @return 
memLinux <- function() {
  temp <- as.numeric(system(paste("ps -p ", Sys.getpid(), " -o rss", sep = ""), intern = TRUE)[2])
  return(temp / 1024)
}

#' @name workers
#' @description 计算合适的并行运算运行核数
#' @param multiCores 逻辑值，是否进行多核运算
#' @param autoCores  逻辑值，是否程序根据预留内存自动分配内核，如果为FALSE，将使用所有可用的线程核
#' @param keepMem    预留百分之多少的内存，防止程序内存不够
#' @return 返回可以使用的核数
workers <- function(multiCores = TRUE, autoCores = TRUE, keepMem = 0.1) {
  library(future)
  gc()
  if (multiCores) {
    if (autoCores) {
      t1 <- memLinux() ## R当前占用内存
      t2 <- sysmem() ## 系统总内存数
      works <- floor(t2 * (1 - keepMem) / t1) ## 可以使用的核数
      if (works < future::availableCores()) {
        works <- works
      } else {
        works <- future::availableCores()
      }
    } else {
      works <- future::availableCores()
    }
  } else {
    works <- 1
  }
  return(works)
}

#' @name logRank
#' @description 计算logrank test p值
#' @name logRank
#' @param data : data.frame格式，临床数据：至少包括Patient_ID,生存时间time,生存状态event,signature分组样本标签sample.label
#' @param time ：数值型，生存时间survival time
#' @param event ：数值型，生存状态survival event
#' @param variate ：因子类型，样本分类标签，sample.label
#' @return return logrank test p value
logRank <- function(data, time, event, variate) {
  options(stringAsFactors = FALSE)
  library(survival)
  
  survObject <- survival::survdiff(survival::Surv(time = time, event = event) ~ variate, data = data)
  p.value <- pchisq(survObject$chisq, length(survObject$n) - 1, lower.tail = F)
  
  return(p.value)
}

#' @name uniVariatePvalue
#' @description 计算单因素cox p值和beta系数
#' @name uniVariatePvalue
#' @param data : data.frame格式，临床数据：至少包括Patient_ID,生存时间time,生存状态event,signature分组样本标签sample.label
#' @param time ：数值型，生存时间survival time
#' @param event ：数值型，生存状态survival event
#' @param variate: 当使用基因表达值作为变量时，则为数值型；当使用样本分类标签为变量时，则为因子型
#' @return 返回一个向量，包含单因素cox beta系数（coef），p值（Pr(>|z|)），HR值（exp(coef)），HR的95%置信区间（lower .95，upper .95） 
uniVariatePvalue <- function(data, time, event, variate) {
  options(stringAsFactors = FALSE)
  library(survival)
  
  # 计算单因素cox p值
  univ.formula <- survival::Surv(time = time, event = event) ~ variate
  univ.model <- survival::coxph(univ.formula, data = data)
  
  Ps <- summary(univ.model)$coefficients[, c("coef", "Pr(>|z|)")]
  HRs <- summary(univ.model)$conf.int[, c("exp(coef)", "lower .95", "upper .95")]
  result <- c(Ps, HRs)
  
  return(result)
}


#' @name multiVariatePvalue
#' @description 计算multivariate cox p值
#' @param data clinical data
#' @param time survival time
#' @param event survival event
#' @param variate character vector, the names of clinical variates
#' @return 返回一个向量，多因素cox p值（Pr(>|z|)），HR值（exp(coef)），HR的95%置信区间（lower .95，upper .95）【只针对于第一个协变量】
multiVariatePvalue <- function(data, time, event, variate) {
  options(stringsAsFactors = FALSE)
  library(survival)
  
  # 计算multivariate cox p值
  multiv.formula <- as.formula(paste("survival::Surv(time = time, event = event) ~", paste(variate, collapse = " + ")))
  multiv.model <- survival::coxph(multiv.formula, data = data)
  
  ## 这里一块儿取出p值和HR，因为单独取p值是没有名字的
  p.value <- summary(multiv.model)$coefficients[grep(pattern = variate[1], rownames(summary(multiv.model)$coefficients)), c("Pr(>|z|)", "exp(coef)")]
  HRs <- summary(multiv.model)$conf.int[grep(pattern = variate[1], rownames(summary(multiv.model)$coefficients)), c("lower .95", "upper .95")]
  result <- c(p.value, HRs)
  
  return(result)
}

#' @name BetaCoefficient
#' @description 计算所有基因的beta系数
#' @param clinical clinical data，至少包含Patient_ID,time,event
#' @param exp.data expression data, 表达矩阵，行为gene，列为样本
#' @param signatures siganture矩阵，每行代表一个signature组合，包含组合中的所有基因
#' @return 返回一个beta系数值的向量，每个元素的名字为系数对应的基因
BetaCoefficient <- function(clinical, exp.data, signatures) {
  options(stringsAsFactors = FALSE)
  
  # 取出signature矩阵中的所有唯一基因
  all.genes <- unique(as.vector(signatures))
  
  # 根据每个基因的表达水平，计算每个基因对应的单因素cox beta系数
  beta <- do.call(c, lapply(all.genes, function(x) {
    uniVariatePvalue(data = clinical, time = clinical$time, event = clinical$event, variate = exp.data[x,])["coef"]
  }))
  names(beta) <- all.genes
  
  return(beta)
}

#' @name RiskScore
#' @description 批量计算每个signature的risk score 
#' @param exp.data expression data, 表达矩阵，行为gene，列为样本
#' @param beta 数值向量，存储的是单因素cox beta值，名字是基因名字,且与exp.data的行名对应
#' @param signatures siganture矩阵，每行代表一个signature组合，包含组合中的所有基因 
#' @return 返回一个risk score矩阵, 每行代表一个signature，每列代表一个样本
RiskScore <- function(exp.data, beta, signatures) {
  
  # 提取signature涉及的所有基因表达值
  all.genes <- unique(as.vector(signatures))
  exp.data <- exp.data[all.genes,]
  
  # 将beta中基因对应的系数值进行排秩，与提取的exp.data中基因顺序一致
  orders <- match(rownames(exp.data), names(beta))
  ordered.beta <- beta[orders]
  
  # 计算基因加权的表达值
  weighted.exp <- exp.data * ordered.beta
  
  # 计算每个signature的risk score得分
  signatures.list <- split(signatures, row(signatures))
  risk.score <- do.call(rbind, furrr::future_map(signatures.list, ~ apply(weighted.exp[.x, ], 2, sum)))
  rownames(risk.score) <- lapply(signatures.list, paste, collapse = "_")
  return(risk.score)
}

#' @name LableSample
#' @description 根据指定的值和和risk score将样本分为高低风险组
#' @param signatrue.risk.score 数值型向量，一个signature的risk score值
#' @param value 数值型，区分高低风险的cutoff值，如果为NULL，自动按risk score中值分类；如果设置为某个数值，则按该数值分高低风险组
#' @return 返回一个向量，标记了样本的高低风险
LableSample <- function(signature.risk.score, value = NULL) {
  
  # 根据signature risk score中值分类样本
  if (is.null(value)) {
    label.sample <- ifelse(signature.risk.score > median(signature.risk.score), "High_risk", "Low_risk")
  } else {
    label.sample <- ifelse(signature.risk.score > value, "High_risk", "Low_risk")
  }
  
  return(label.sample)
}


