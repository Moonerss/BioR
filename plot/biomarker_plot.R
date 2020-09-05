#------------------------------------------------------------
#                  |Expression Biomarker|
#
#           Author: Erjie Zhao
#           Time :  2020-08-13
#          Update:
#           绘制表达riskscore高低风险的三联图
#------------------------------------------------------------

#' @param dat 输入的数据必须包含patient, time,event,剩余列为基因表达数据或其他数据列
#' @param cutoff 区分高低风险的方式
#' @param genes 模型中包含的基因（必须是dat的列名）
#' @param weight 模型中基因的权重（必须与genes长度一致）
#' @import survival, ggplot2, ggplotify, gridExtra

riskscore_plot <- function(dat, cutoff = c("median", "mean"), genes = NULL, weight = NULL) {

  ## 不提供特定基因，直接计算剩余列的所有基因  
  if (is.null(genes)) {
    gene_name <- setdiff(colnames(dat), c("patient", "time", "event"))
  } else {
    gene_name <- genes
  }

  ## 构建Surv对象
  if (length(gene_name) == 1) {
    if (is.null(weight)) {
      surv_obj <- eval(prase(paste("Surv(time, event) ~", gene_name)))
    } else {
      surv_obj <- eval(prase(paste("Surv(time, event) ~", weight , "*", gene_name)))
    }
  } else {
    if (is.null(weight)) {
      surv_obj <- eval(prase(paste("Surv(time, event) ~", paste(gene_name, collapse = " + "))))
    } else {
      surv_obj <- eval(prase(paste("Surv(time, event) ~", paste(paste(weight, gene_name, sep = "*"), collapse = " + "))))
    }
  }

  ## 构建cox模型
  model <- coxph(surv_obj, data = dat)

  ## 使用Predict函数，计算出每位患者的risk score
  RiskScore <- predict(model, type = "risk")
  names(RiskScore) = dat$patient

  ## 开始绘制风险模型的risk score 点图
  fp <- RiskScore
  phe <- dat
  fp_dat <- data.frame(patientid = 1:length(fp), fp = as.numeric(sort(fp)))

  ## 添加风险分组，以风险评分的中位值将患者分为两组
  if (cutoff == "median") {
    fp_dat$riskgroup <- ifelse(fp_dat$fp >= median(fp_dat$fp),'high','low')
  } else if (cutoff == "mean"){
    fp_dat$riskgroup <- ifelse(fp_dat$fp >= mean(fp_dat$fp),'high','low')
  }

  ## 绘制生存事件的点图
  sur_dat <- data.frame(patientid = 1:length(fp), time = phe[names(sort(fp)),'time'], event = phe[names(sort(fp )),'event'])
  sur_dat$event <- ifelse(sur_dat$event==0,'alive','death')
  sur_dat$event <- factor(sur_dat$event,levels = c("death","alive")) 

  ## 表达热图
  exp_dat <- dat[names(sort(fp)), genes]

  ## risk score 排布图
  p1 <- ggplot(fp_dat, aes(x = patientid, y = fp)) +
        geom_point(aes(color = riskgroup)) + 
        scale_colour_manual(values = c("red","green")) + 
        theme_bw() +
        labs(x = "Patient ID(increasing risk score)", y = "Risk score") +
        geom_hline(yintercept = median(fp_dat$fp), colour = "black", linetype = "dotted", size = 0.8) +
        geom_vline(xintercept = sum(fp_dat$riskgroup == "low"), colour="black", linetype = "dotted", size = 0.8)
  
  ## 生存时间样本排布图
  p2 <- ggplot(sur_dat, aes(x = patientid, y = time)) +
        geom_point(aes(col = event)) +
        theme_bw() +
        scale_colour_manual(values = c("red","green")) +
        labs(x = "Patient ID(increasing risk score)", y = "Survival time(year)") +
        geom_vline(xintercept = sum(fp_dat$riskgroup=="low"),colour = "black", linetype = "dotted", size = 0.8)

  ## 表达等的聚类热图
  mycolors <- colorRampPalette(c("white", "green", "red"), bias = 1.2)(100)
  tmp <- t(scale(exp_dat))
  tmp[tmp > 1] <- 1
  tmp[tmp < -1] <- -1
  p3 <- pheatmap(tmp,col = mycolors, show_colnames = F,cluster_cols = F)

  ## 拼接图形
  plots <- list(p1,p2,as.ggplot(as.grob(p3)))
  lay1 <- rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7))) # 图片布局
  res <- grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,3,2),weights=c(10,10,10))

  return(res)
}



#------------------------------------------------------------
#                  |R Global Setting|
#
# You can (un)comment any code you dislike.
#   Any Question, please email to
#       Shixiang Wang <w_shixiang@163.com>
#   or file an issue to
#       <https://github.com/ShixiangWang/MessageBoard/issues>
#------------------------------------------------------------


# Global options ----------------------------------------------------------
options(encoding = "UTF-8") # Set file encoding
Sys.setenv("LANGUAGE" = "EN") # Set language of R message

# Package download mirrors ------------------------------------------------
## For Bioconductor packages
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
## For CRAN packages
## Full list see mirrors on <https://cran.r-project.org/>
options("repos" = c(CRAN = "https://mirrors.tongji.edu.cn/CRAN/"))

# Set local custom R library path -----------------------------------------
.CUSTOM_LIB <- "~/Library/R" # Set your custom library location
# Please do not add '/' at the end !!!

if (!dir.exists(.CUSTOM_LIB)) {
  dir.create(.CUSTOM_LIB, recursive = TRUE)
}

.libPaths(c(.CUSTOM_LIB, .libPaths()))
message("Using library: ", .libPaths()[1])

# Set R temp directory ----------------------------------------------------
## Uncomment the following code if you want to set R temp directory
# .TMP = "~/Rtmp"
# if(dirname(tempdir()) != .TMP){
#   if(!dir.exists(.TMP)) dir.create(.TMP)
#   cat(paste0("TMPDIR = ", .TMP), file="~/.Renviron", sep = "\n")
# }
# message("Using temp directory: ", .TMP)


# Use pacman to manage R packages -----------------------------------------
# Reference: <https://www.jianshu.com/p/cb16ded75672>
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", dependencies = TRUE)
}

library(pacman)

# VSCode plugin setting ---------------------------------------------------
# Only use when you code R with VSCode
# Reference: <https://github.com/REditorSupport/vscode-r-lsp>
if (!requireNamespace("languageserver", quietly = TRUE)) {
  pacman::p_install(languageserver)
}

source(file.path(
  if (.Platform$OS.type == "windows") {
    file.path(Sys.getenv("HOMEDRIVE"), Sys.getenv("HOMEPATH"))
  } else {
    Sys.getenv("HOME")
  }, ".vscode-R", "init.R"
))

# Global utilities --------------------------------------------------------

## Try installing R packages again and again
## It may be useful when install GitHub R packages
.loop_install <- function(pkgs, force = FALSE) {
  ## If force=TRUE,
  ## force installation when packages already exist on local system.
  Sys.sleep(1)
  tryCatch(
    {
      message("=> Try installing ", paste(pkgs, collapse = ", "))
      gh_pkg <- pkgs[grepl("/", pkgs)]
      ot_pkg <- setdiff(pkgs, gh_pkg)
      if (length(ot_pkg) > 0) {
        pacman::p_install(ot_pkg, character.only = TRUE, force = force)
      } else {
        if (!requireNamespace("remotes", quietly = TRUE)) {
          pacman::p_install(remotes)
        }
        remotes::install_github(gh_pkg, force = force)
      }
    },
    error = function(e) {
      .loop_install(pkgs, force)
    }
  )
}

## Utilities from rvcheck, xfun and other R package
## Thanks to the authors
if (!requireNamespace("rvcheck", quietly = TRUE)) {
  pacman::p_install(rvcheck)
}
if (!requireNamespace("xfun", quietly = TRUE)) {
  pacman::p_install(xfun)
}

## Load function from package
.get_fun <- rvcheck::get_fun_from_pkg
## Open any directory in any operating system
.open <- rvcheck::o
## Check whether packages were installed
.is_installed <- rvcheck::is.installed
## Check if on a R server
.is_rserver <- suppressMessages(.get_fun("rvcheck", "is.rserver"))
## Simpler function to download file
.download_file <- xfun::download_file
.upload_ftp <- xfun::upload_ftp
## Functions to obtain, remove, and change extensions in filenames
.fl_ext <- xfun::file_ext
.rm_ext <- xfun::sans_ext
.ch_ext <- xfun::with_ext
## Provide the "file" version of gsub(),
## i.e., they perform searching and replacement in files via gsub().
.gsub_file <- xfun::gsub_file
.gsub_files <- xfun::gsub_files
.gsub_dir <- xfun::gsub_dir
.gsub_ext <- xfun::gsub_ext
## Change the working directory, evaluate the expression,
## and restore the working directory
.move_run <- xfun::in_dir
## Install a source package from a directory
.install_dir <- xfun::install_dir
## Check OS
.is_linux <- xfun::is_linux
.is_windows <- xfun::is_windows
.is_macos <- xfun::is_macos
.is_unix <- xfun::is_unix
## Number to words, e.g. 1 to one, 10 to ten
.numbers_to_words <- xfun::numbers_to_words