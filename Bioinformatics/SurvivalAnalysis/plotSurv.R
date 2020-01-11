##  TO DO: 绘制生存曲线图
##
##  Author: jeason zhao
##################################################


#' @name plotSurv
#' @description 本函数基于survminer包完成一系列标注  
#' @param clinical.data 临床数据，至少包含四列：Patient_ID, event, time, sample.label（factor类型）
#' @param label.order 字符向量，为防止sample.label未设置成factor，可以在此设置factor排序，若sample.label为factor，则不进行转换。默认为NULL
#' @param upper.time 数值型，生存时间上限（与time单位保持一致），默认为NULL，若使用该参数，则删除超过年限的样本
#' @param xscale 字符，负责x轴时间单位转换。允许的选项包括"d_m"，"d_y"，"m_d"，"m_y"，"y_d"和"y_m"。xscale =“d_m”会将x轴单位从天转换为月
#' @param xlab 字符，x轴的标签
#' @param median.time 逻辑值，是否绘制中位生存时间和95%置信区间
#' @param surv.median.line 中位生存时间的标记方式，允许的值包括c（"none"，"hv"，"h"，"v"）
#' @param HR 逻辑值，是否在图中标注出HR信息。注意：当且仅当样本分为两组时可以使用
#' @param risk.table 逻辑向量，是否绘制risk.table，默认为TRUE
#' @param conf.int 逻辑向量，是否画出置信区间，默认为FALSE
#' @param pval 逻辑向量，是否在图中标出log rank p值，默认为TRUE
#' @param ylab 字符，y轴的标签
#' @param main 字符，主标题的名字
#' @return 返回一个ggplot对象
#' 

plot.surv <- function(clinical.data, label.order = NULL, upper.time = NULL,
                      xscale = 1, xlab = "Time", median.time = TRUE, 
                      surv.median.line = "none", HR = FALSE, risk.table = TRUE,
                      pval = TRUE, conf.int = FALSE, main = NULL,
                      ylab = "Survival probability") {
  
  #载入相关R包
  require(survival)
  require(survminer)
  require(RColorBrewer)
  require(gridExtra)
  
  
  #如果设置upper.time，则去除生存时间超过upper.time的样本
  if (!is.null(upper.time)) {
    clinical.data <- clinical.data[clinical.data$time <= upper.time,]
  }
  
  # 判断样本标签的类型，集体转换成factor类型
  if (!is.factor(clinical.data$sample.label)) {
    if (is.null(label.order)) {
      cat("sample.label 不是factor类型，程序自动进行转换...\n")
      clinical.data$sample.label <- as.factor(clinical.data$sample.label)
    } else {
      clinical.data$sample.label <- factor(clinical.data$sample.label, levels = label.order)
    }   
  } else {
    cat("sample.lable 是factor类型，不进行转换...\n")
  }
  
  t.name <- levels(clinical.data$sample.label)
  
  cat("factor levels :", t.name, "\n")
  
  ## 分配默认颜色
  if (length(t.name) > 6) {
    stop("样本分组>6，超过函数接受范围")
  }
  colors <- c("#808080","#EA4335","#4285F4","#FBBC05","#34A853","#000000") # 顺序：灰，红，蓝，黄，绿，黑
  t.col <- colors[seq(length(t.name))]
  
  # 构造生存对象
  km.curves <- survfit(Surv(time, event) ~ sample.label, data = clinical.data)
  
  # 计算HR值和95%CI
  if (length(t.name) == 2) {
    if (HR) {
      cat("All compared with the ", t.name[1], "in cox analysis...\n")
      cox.obj <- coxph(Surv(time, event) ~ sample.label, data = clinical.data)
      tmp <- summary(cox.obj)
      HRs <- round(tmp$coefficients[ ,2], digits = 2)
      HR.confint.lower <- round(tmp$conf.int[,"lower .95"], 2)
      HR.confint.upper <- round(tmp$conf.int[,"upper .95"], 2)
      HRs <- paste0(HRs, " (", HR.confint.lower, "-", HR.confint.upper, ")")	  
    }
  }
  
  # 构造生存图像中图例显示文字
  legend.content <- substr(names(km.curves$strata), start = 14, stop = 1000)
  
  # x轴刻度单位转换
  if (is.numeric(xscale) | (xscale %in% c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"))) {
    xscale = xscale
  } else {
    stop('xscale should be numeric or one of c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m").')
  }
  
  # 隐函数：转换生存时间单位
  .format_xticklabels <- function(labels, xscale){
    # 1 year = 365.25 days
    # 1 month = 365.25/12 = 30.4375 days
    if (is.numeric(xscale)) {
      xtrans <- 1/xscale
    } else {
      xtrans <- switch(xscale,
                       d_m = 12/365.25,
                       d_y = 1/365.25,
                       m_d = 365.25/12,
                       m_y = 1/12,
                       y_d = 365.25,
                       y_m = 12,
                       1
      )
      round(labels*xtrans,2)	
    }
  }
  
  # 在图中添加中位生存时间及其95%CI,放在副标题位置
  subtitle <- NULL
  if (median.time) {
    if (is.numeric(xscale)) {
      median.km.obj <- km.curves
    } else if (xscale %in% c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m")) {
      clinical.data$time <- .format_xticklabels(labels = clinical.data$time, xscale = xscale)
      median.km.obj <- survfit(Surv(time, event) ~ sample.label, data = clinical.data)
    }
    survival.time.info <- NULL
    survival.time.info <- rbind(survival.time.info, summary(median.km.obj)$table)
    median.survival <- round(survival.time.info[!duplicated(survival.time.info[,7:9]),7:9], digits = 2) # 注意：这里取得的置信区间上界可能为NA
    if (length(levels(clinical.data$sample.label)) == 1) {
      tmp1 <- levels(clinical.data$sample.label)
    } else {
      tmp1 <- do.call(rbind,strsplit(rownames(summary(median.km.obj)$table), split = "="))[,2]
    }
    tmp2 <- paste(median.survival[,1], "(", median.survival[,2], "-", median.survival[,3], ")")
    subtitle <- paste(tmp1, tmp2, sep = ":", collapse = "\n")
  }
  
  # ggsurvplot绘制生存图像
  ggsurv <- ggsurvplot(km.curves,               # survfit object with calculated statistics.
                       data = clinical.data,             # data used to fit survival curves.
                       palette = t.col,
                       
                       #图的主题构架
                       risk.table = risk.table,       # show risk table.
                       pval = pval,             # show p-value of log-rank test.
                       surv.median.line = surv.median.line,  # add the median survival pointer.
                       title = main,     #主标题名字
                       subtitle = subtitle, #副标题
                       font.main = 15,       #主标题字体大小              
                       xlab = xlab,   # customize X axis label.
                       ylab = ylab,   # customize Y axis label
                       xscale = xscale,
                       
                       
                       #图例设置
                       legend.title = "", #图例标题，一般不用，设为空
                       legend.labs = legend.content, #图例文字描述
                       legend = c(0.8,0.9), #图例的位置，取值在【0,1】之间
                       font.legend = 9,     #图例字体大小
                       
                       #risk table设置
                       tables.theme = theme_cleantable(),#table主题
                       risk.table.title = "No. at risk:",#table标题
                       risk.table.y.text.col = T, # 使用颜色代替Y轴文字
                       risk.table.y.text = FALSE, # Y轴不使用文字注释
                       tables.height = 0.15,      # table的高度
                       risk.table.fontsize = 3    # risk table内文字的大小
  )
  # HR标注的位置限定
  if (length(t.name) == 2) {
    if (HR) {
      ggsurv$plot <- ggsurv$plot + ggplot2::annotate("text", x = max(km.curves$time)/13, 
                                                     y = 0.15, size = 5, label = paste("HR=", HRs))	  
    }
  } 
  #图的标题居中
  ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(size = 10), 
                                     plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))
  #表的标题
  ggsurv$table <- ggsurv$table + theme(plot.title = element_text(hjust = -0.04),
                                       plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))
  
  # 判断分类的类数，如果只有两类，就不必计算两两之间的log rank p值
  if(length(t.name) > 2) {
    # 计算pairwise的log rank的p值
    res <- pairwise_survdiff(Surv(time, event)~sample.label, data=clinical.data);
    pairwise.pvalue <- round(res$p.value, digits = 4);
    pairwise.pvalue[which(pairwise.pvalue < 0.0001)] <- "<0.0001"
    pairwise.pvalue[is.na(pairwise.pvalue)] <- "-"
    
    # 添加表格
    tt <- ttheme_minimal(core = list(fg_params = list(col = "black"),bg_params = list(fill = NA, col = "black")),
                         colhead = list(fg_params = list(col = NA),bg_params = list(fill = t.col, col = "black")),
                         rowhead = list(fg_params = list(col = NA, hjust = 1),bg_params = list(fill = c("white",t.col[-1]), col = "black")))
    pairwise.table <- tableGrob(pairwise.pvalue, theme = tt)
    ggsurv <- ggarrange(ggarrange(ggsurv$plot, ggsurv$table, nrow=2, heights=c(2,0.5)),
                        pairwise.table, nrow=2, heights = c(2,0.5),
                        labels = c("","p from pairwise comparisons"),
                        hjust = 0, font.label = list(size = 15, face = "plain"))
  }
  
  ggsurv
}
