#' @@ 预测子相关性热图
#' @param  scoreMatrix 预测子score矩阵，行为样本，列为预测子
#' @param  method 指定可视化的方法，可以是圆形、方形、椭圆形、数值、阴影、颜色或饼图形
#' @param  type 指定展示的方式，可以是完全的、下三角或上三角
#' @param  addrect 当order为hclust时，可以为添加相关系数图添加矩形框，默认NULL(不添加框)，如果想添加框时，只需为该参数指定一个整数即可
#' @param  tl.cex 文本标签(变量名称)的大小，默认的 0.75倍
#' @param  number.cex 写在图中的相关系数文字的大小，默认的0.7倍
#' @param  adjust 计算相关性时，是否进行多重检验矫正，默认方法有 "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param  title 题目
#' @returnType 
#' @return 
#' 
#' @author zhaoej

        

correlation.heatmap <- function(scoreMatrix,
                                method = c("circle", "square", "ellipse", "number", "shade", "color", "pie"), 
                                type = c("full", "lower", "upper"), 
                                addrect = NULL, tl.cex = 0.6, number.cex = 0.4, adjust = c("fdr", "none","bonferroni", "BH"),
                                title = "",
                                addCoef.col = "#666666", addgrid.col = NA){
  
        method <- match.arg(method)
        type <- match.arg(type)
        adjust <- match.arg(adjust)
        
        library(psych)
        res <- corr.test(scoreMatrix, method = "spearman", adjust = adjust)
        # 相关性p值，对角线下方为原始p值，上方为矫正后的p值
        a <- res$p
        a[lower.tri(a)] <- 0
        b <- t(a)
        p.matrix <- a+b
        diag(p.matrix) <- 1
        
        library(corrplot)
        col3 <- colorRampPalette(c("#0066CC", "white", "#FF6666")) 
        
        # diag：是否展示对角线上的结果，默认为TRUE
        # addgrid.col: grid的颜色，如果等于NA，不添加grid
        # addCoef.col: 图中相关系数的颜色，如果等于NULL，不写相关系数
        corrplot(res$r, method = method, type = type, order = "hclust", addrect = addrect, 
                 addgrid.col = addgrid.col, addCoef.col = addCoef.col,
                 col = col3(20), p.mat = p.matrix, insig = "label_sig", sig.level = 0.05, pch.cex = 2, pch.col = "#CCCCCC",  tl.cex = tl.cex, rect.col = "black", number.cex = number.cex, tl.col = "black", diag = F,
                 title = title)
}
