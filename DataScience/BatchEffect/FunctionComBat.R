#ComBat是基于经典贝叶斯的分析方法，运用已知的批次信息对高通量数据进行批次校正。
#在sva R package 中提供了ComBat用于处理批次效应。ComBat有两个方法可供选择，一种是基于参数和一种非参数方法，combat函数的par.prior参数可以设置。
#函数输入数据为经过标准化的数据矩阵，返回结果为经过批次校正后的一个数据矩阵
# 参考文献：Adjusting batch effects in microarray expression data using empirical Bayes methods. 2007 Biostatistics

#' @@ 矫正数据的批次效应
#' @param expressionMatrix matrix格式表达谱，不能为data.frame，The input data are assumed to be cleaned and normalized before batch effect removal，行为基因，列为样本。
#' 一般最好输入为log后表达矩阵，要求应将在样本中全为0的基因删除，否则可能会产生NaN 
#' @param batch 数值型向量，指明批次，要求和表达谱中样本名一一对应，如1,1,1,2,2,2。
#' @param covariate 因子型data.frame,每一列为一个covariates(感兴趣变量),要求为因子型,且和batch对应的样本顺序一致. 如去除不同实验平台的批次,计算正常和癌症的差异基因,则每个实验中的normal和cancer样本分组则为感兴趣的变量,即covariate = c("normal", "cancer", "normal", "normal",...)
# 协变量是一般为后续研究感兴趣的变量,批次矫正时需要考虑在各个批次中数量分布情况差异(不平衡情况)
#' @param par.prior TRUE/FALSE，TRUE表示使用parametric adjustments方式，FALSE表明使用non-parametric adjustments方式。By default, it performs parametric empirical Bayesian adjustments. 
#' If you would like to use nonparametric empirical Bayesian adjustments, use the par.prior=FALSE option (this will take longer)
#' @param prior.plots TRUE/FALSE ,如果为TURE，则会绘制EB批次效应参数γig和δig2的 kernel density estimate of the empirical values（黑线）和EB-based prior distribution（红线）

#' @returnType matrix or data.frame
#' @return 矫正批次后的表达矩阵
#' 
#' @author 龙志林 2019.1.7

#' @note 
#使用combat矫正批次会产生负数，sva package的作者建议输入的表达矩阵为log后的标化矩阵，并且可以将输出的矩阵小于0的值设为0 --- 参考https://support.bioconductor.org/p/73266/；https://www.biostars.org/p/173684/；

# 在sva包中，假定有两种变量需要考虑：1.兴趣变量（covariates, 如癌症和正常对照）。2.调整变量,一般指批次
# 另外有两种模型矩阵（model matrices):1. full model(全模型):包含以上的两种变量；2. null model：只包含调整变量，即需要矫正的变量。
# 基于sva作者Johnson建议(https://www.biostars.org/p/196430/):  
# covariates should ALWAYS be added unless you have a completely balanced factorial design, in which it won't make much difference whether you add them or not. However, if you have an unbalanced proportion of treatment/controls or if some treatments are missing from one of the batches, you MUST include covariates or you risk introducing bias or removing biological signal from the data.

# @example
# library(sva)
# library(bladderbatch)
# data(bladderdata)
# pheno = pData(bladderEset)
# expressionMatrix = exprs(bladderEset)
# batch = pheno$batch
# covariate = data.frame(cancer = as.factor(pheno$cancer))
# correctMatrix <- FunctionCombat(expressionMatrix, batch, covariate = covariate)

FunctionCombat <- function(expressionMatrix, batch, covariate = NULL, par.prior = TRUE, prior.plots = FALSE){
        library(sva)
        # if(is.factor(batch)){
        #         batch <- batch
        # }else{
        #         batch <- as.factor(batch)  
        # }
        
        # 由于当前只进行batch的矫正，所以建立null model
        phen <- data.frame(sample = 1:ncol(expressionMatrix), batch = batch)
        rownames(phen) <- colnames(expressionMatrix)
        
        if(is.null(covariate)){
                # 感兴趣变量缺失或不清楚,这种模式可能会丢失部分有意义的关系
                mod <- model.matrix(~1, data = phen)
        }else{
                # 非平衡实验设计
                dataType <- sapply(covariate, as.factor)
                # if(length(grep("FALSE", dataType)) > 0){
                #         stop("Covariate parameter is non-factor types!")
                # }
                mod <- model.matrix(~., data = covariate)
        }
        
        # 矫正数据批次
        combat_edata <- ComBat(dat=expressionMatrix, batch=batch, mod=mod, par.prior=par.prior, prior.plots = prior.plots)
        
        # 将矫正后的表达值为负数的设置为0
        combat_edata[which(combat_edata<0)] <- 0
        
        return(combat_edata)
}
