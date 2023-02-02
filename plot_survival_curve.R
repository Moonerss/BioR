#######################################################################
## Filename: plot_survival_curve.R
## Created: 2021-10-11
##
## Purpose: plot survival curve with many other element
##
#######################################################################


#######################################################################################################
##
##                                            plot survival curve
## Argument:
##       clinical_data         - clinical data to do survival analysis, at least three columns: 
##                               time, event, group
##       upper_time            _ a numeric value as Upper limit of survival time, which must be 
##                               consistent with the unit in clinical.data, NULL by default.
##                               If used, samples with survival time greater than upper_time 
##                               will be removed.
##       xscale                _ change the unit of x axis, can use "d_m", "d_y", "m_d", "m_y",
##                               "y_d", and "y_m". y = 'year', m = 'month', d = 'day'
##       main                  - plot main title
##       xlab                  - x axis title
##       ylab                  - y lab title
##       median_time           - whether show median time in subtitle
##       median_line           - character vector for drawing a horizontal/vertical line at median
##                               survival. Allowed values include one of c("none", "hv", "h", "v").
##                               v: vertical, h:horizontal. 
##       show_risk_table       - whether show risk table 
##       show_pval             - whether show p value
##       show_HR               - whether show HR
##       conf_int              - whether show confident interval in subtitle
##       show_subtitle         - whether show subtitle
##       show_pairwise_table   - whether show pairwise log rank p value table
##
#######################################################################################################
plot_survival_curve <- function(clinical_data, upper_time = NULL, xscale = NULL,
                                main = NULL, xlab = 'Time', ylab = 'Survival probability',
                                median_time = TRUE, median_line = 'none',
                                show_risk_table = TRUE, show_pval = TRUE, show_HR = FALSE,
                                conf_int = FALSE, show_subtitle = FALSE, show_pairwise_table = FALSE) {
  # attach packages
  require(survival)
  require(survminer)
  require(RColorBrewer)
  require(gridExtra)
  require(dplyr)
  require(magrittr)
  
  # if set upper_time, remove samples with survival time > upper_time
  if (!is.null(upper_time)) {
    clinical_data <- clinical_data %>% filter(time <= upper_time)
  }
  
  ## constrict group color
  if (!is.factor(clinical_data$group)) {
    clinical_data$group <- as.factor(clinical_data$group)
  }
  
  group_name <- levels(clinical_data$group)
  
  if (length(group_name) > 6) stop("More than 6 groups exist!")
  colors <- c("#808080","#EA4335","#4285F4","#FBBC05","#34A853","#000000") # 顺序：灰，红，蓝，黄，绿，黑
  group_col <- colors[1:length(group_name)]
  
  # survfit
  km_curves <- survfit(Surv(time, event)~group, data = clinical_data)
  
  # calculate HR and 95%CI
  if (length(group_name) == 2) {
    if (show_HR) {
      cox.obj <- coxph(Surv(time, event)~group, data = clinical_data)
      tmp <- summary(cox.obj)
      HRs <- round(tmp$coefficients[ ,2], digits = 2)
      HR.confint.lower <- round(tmp$conf.int[,"lower .95"], 2)
      HR.confint.upper <- round(tmp$conf.int[,"upper .95"], 2)
      HRs <- paste0(HRs, " (", HR.confint.lower, "-", HR.confint.upper, ")")	  
    }
  }
  
  # content legend
  legend.content <- substr(names(km_curves$strata), start = 14, stop = 1000)
  
  # the scale of x axis
  if (is.null(xscale)) {
    xscale <- 1
  }
  if (is.numeric(xscale) | (xscale %in% c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"))) {
    xscale = xscale
  } else {
    stop('xscale should be one of c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m").')
  }
  
  # transform time unit
  .format_xticklabels <- function(labels, xscale){
    # 1 year = 365.25 days
    # 1 month = 365.25/12 = 30.4375 days
    if (is.numeric(xscale)) xtrans <- 1/xscale
    else
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
  
  ## add median survival time and 95% CI in subtitle
  subtitle <- NULL
  if (show_subtitle) {
    if (median_time) {
      if (is.numeric(xscale)) {
        median.km.obj = km_curves
      } else if (xscale %in% c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m")) {
        clinical.data$time <- .format_xticklabels(labels = clinical_data$time, xscale = xscale)
        median.km.obj <- survfit(Surv(time, event)~group, data = clinical_data)
      }
      survival.time.info <- NULL
      survival.time.info <- rbind(survival.time.info, summary(median.km.obj)$table)
      median.survival <- round(survival.time.info[!duplicated(survival.time.info[,7:9]),7:9], digits = 2) # 注意：这里取得的置信区间上界可能为NA
      if (length(levels(clinical_data$group)) == 1) {
        tmp1 <- levels(clinical_data$group)
      } else {
        tmp1 <- do.call(rbind,strsplit(rownames(summary(median.km.obj)$table), split = "="))[,2]
      }
      tmp2 <- paste(median.survival[,1], "(", median.survival[,2], "-", median.survival[,3], ")")
      subtitle <- paste(tmp1, tmp2, sep = ":", collapse = "\n")
    }  
  }
  
  ## ggsurvplot function
  ggsurv <- survminer::ggsurvplot(km_curves,                # survfit object with calculated statistics.
                                  data = clinical_data,     # data used to fit survival curves
                                  palette = group_col,
                                  # main of the theme
                                  risk.table = show_risk_table,   # show risk table.
                                  pval = show_pval,               # show p-value of log-rank test.
                                  surv.median.line = median_line, # add the median survival pointer.
                                  title = main,
                                  subtitle = subtitle,
                                  font.main = 15,
                                  xlab = xlab,                    # customize X axis label.
                                  ylab = ylab,                    # customize Y axis label.
                                  xscale = xscale, 
                                  # legend
                                  legend.title = '',
                                  legend.title = legend.content,
                                  legend = c(0.8, 0.9),           # position of legend.
                                  font.legend = 9,
                                  # risk table
                                  tables.theme = theme_cleantable(), # table theme.
                                  risk.table.title = "No. at risk:", # table title.
                                  risk.table.y.text.col = T,         # use color replace text.
                                  risk.table.y.text = FALSE, 
                                  tables.height = 0.15,      # table height.
                                  risk.table.fontsize = 3    # risk table text size.
                                  )
  # HR height position
  if (length(group_name) == 2) {
    if (show_HR) 
      ggsurv$plot <- ggsurv$plot + ggplot2::annotate("text", x = max(km_curves$time)/13, 
                                                     y = 0.15, size = 5, label = paste("HR=", HRs))
  } 
  
  # plot title position
  ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(size = 10), 
                                     plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))
  # table title
  ggsurv$table <- ggsurv$table + theme(plot.title = element_text(hjust = -0.04),
                                       plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))
  
  # show pairwise log rank p
  if (show_pairwise) {
    if(length(group_name) > 2) {
      # pairwise log rank p value
      res <- survminer::pairwise_survdiff(Surv(time, event)~group, data = clinical_data);
      pairwise.pvalue <- round(res$p.value, digits = 4);
      pairwise.pvalue[which(pairwise.pvalue < 0.0001)] <- "<0.0001";
      pairwise.pvalue[is.na(pairwise.pvalue)] <- "-"
      
      # add p value table
      tt <- ttheme_minimal(core = list(fg_params = list(col = "black"),bg_params = list(fill = NA, col = "black")),
                           colhead = list(fg_params = list(col = NA),bg_params = list(fill = group_col, col = "black")),
                           rowhead = list(fg_params = list(col = NA, hjust = 1),bg_params = list(fill = c("white", group_col[-1]), col = "black"))
      )
      pairwise.table <- tableGrob(pairwise.pvalue, theme = tt)
      ggsurv <- ggarrange(ggarrange(ggsurv$plot, ggsurv$table, nrow=2, heights=c(2,0.5)),
                          pairwise.table, nrow=2, heights = c(2,0.5),
                          labels = c("","p from pairwise comparisons"),
                          hjust = 0, font.label = list(size = 15, face = "plain"))
      
    }
  } else {
    ggsurv <- ggarrange(ggsurv$plot, ggsurv$table, nrow=2, heights=c(2,0.5))
  }
  
  return(ggsurv)
}


