### 一元线性回归分析函数：univLm
#' @description 对数据框内多列特征分别进行一元线性回归分析，以表格形式输出结果
#' @param dependent 响应变量列名，数值型向量
#' @param control 控制变量，默认为NULL
#' @param data 行为样本，列为特征和响应变量的数据框
#' @return 返回一个数据框，包含一元线性回归分析的结果信息, 信息包括单多因素的OR（95% CI for HR），p value
univLm <- function(data, dependent, control = NULL)
  {
  ###设置工作环境
  options(stringsAsFactors = FALSE, warn = -1);
  
  ###判断自变量类型：数值型（num_covariate），非数值型(nonnum_covariate)，便于后续输出
  covariates <- setdiff(colnames(data), append(dependent, control))
  covariates_type <- sapply(data[,covariates], class)
  
  # num.variate <- NULL
  # for(i in covariates)
  # {
  #   if(is.numeric(clinical.data[, i]))
  #   {
  #     num.variate <- append(num.variate, i)
  #   }
  # }
  # chara.variate <- setdiff(covariates, num.variate) #字符型变量
  
    #STEP1:构建单因素分析的对象
    if(is.null(control)){
      univ_formulas <- sapply(covariates, function(x) as.formula(paste0('dependent~', x)))
    }else{
      univ_formulas <- sapply(covariates, function(x) as.formula(paste0('dependent~', paste(control, collapse = "+"), "+", x)))
    }
    
    #STEP2: 一元线性回归分析
    univ_models <- lapply(univ_formulas, function(x){lm(x, data = data)});
    
    #STEP3:提取有用信息
    univ_results_list <- lapply(univ_models, function(x)
    {                             
      tmp <-summary(x);
      
      #提取p值，保留两位有效数字
      p.value <- round(tmp$coefficients[-1,'Pr(>|t|)'], digits = 4); ##其中-1是为了去除常数项
      #p.value[which(p.value < 0.0001)] <- "<0.0001";
      
      #提取beta值，这里的coefficients为矩阵，但是只有一行，所以可以这样取值
      beta <- round(tmp$coefficients[-1,'Estimate'], digits = 4);
      #提取标准误差
      #SE <- round(tmp$coefficients[-1,'Std. Error'], digits = 4);
      #提取t值
      #Tvalue <- round(tmp$coefficients[-1,'t value'], digits = 4);
      #提取风险比
      OR <- round(exp(coef(x))[-1], digits = 4); ##其中-1是为了去除常数项
      
      #提取95%置信区间上下界
      OR.confint.lower <- round(exp(confint(x))[-1,1], 4); ##其中-1是为了去除常数项
      OR.confint.upper <- round(exp(confint(x))[-1,2], 4); ##其中-1是为了去除常数项 
      
      #合并风险比HR和置信区间为一个内容
      OR <- paste0(OR, " (", OR.confint.lower, "-", OR.confint.upper, ")");
      
      variate <- rownames(tmp$coefficients)[-1];
      
      #将所有值合并在一个矩阵中
      all.data <- as.data.frame(cbind(variate, beta, OR, p.value));
    });
    
    #STEP4:标准化输出格式
    #对于数值型自变量
    # for(i in names(covariates_type)[which(covariates_type == "numeric")])
    # {
    #   tmp <- univ_results[[i]];
    #   tmp$type <- " ";
    #   univ_results[[i]] <- tmp; 
    # }
    # 
    # for(i in names(covariates_type)[which(covariates_type != "numeric")])
    # {
    #   tmp <- univ_results[[i]]
    #   tmp.variate <- substr(tmp$variate, start = 0, stop = nchar(i));
    #   tmp.type <- sub(pattern = i, replacement = "", tmp$variate);
    #   tmp$variate <- tmp.variate;
    #   tmp$type <- paste(tmp.type, "VS", setdiff(unique(data[, i]), tmp.type));
    #   univ_results[[i]] <- tmp; 
    # }
    # 
    # #STEP5: list转化为数据框输出 	
    # univ.result <- do.call(rbind.data.frame, univ_results);
    # univ.result <- univ.result[,c(1,4,2,3)];
    # colnames(univ.result) <- c('variate', 'type', 'univ OR (95% CI for HR)', 'univ p value')
    # return(univ.result)  
    
    
    #STEP5: list转化为数据框输出，每个list元素只保留关心的特征（最后一行），作为控制变量的特征结果删掉
    univ_results <- sapply(univ_results_list, function(x){
      x[nrow(x),]
    })
    univ_results <- t(univ_results)
    univ_results <- data.frame(variate = names(univ_models), beta = as.numeric(univ_results[,'beta']),
                               OR = as.character(univ_results[,'OR']), p_value = as.numeric(univ_results[,'p.value']))
    return(univ_results)
  }