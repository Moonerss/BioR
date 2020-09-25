#' @@ PCA分析
#' @param expressionMatrix matrix or data.frame，表达矩阵，行为基因，列为样本 
#' @param groupLabel factor，分组变量，指明数据中的分组信息 
#' @param plotFile 字符变量，pdf的绝对存储地址，如 /home/longzl2017/dataScience/DimensionalityReduction/test.pdf
#' @param is.log TRUE/FALSE，指明表达谱是否log化，默认为TURE
#' @param label 指明图形上点是否给对应出样本名。
#' Allowed values are "none". 
#' "ind" can be used to label only active individuals. 
#' "ind.sup" is for supplementary individuals. 
#' "quali" is for supplementary qualitative variables. 
#' "var" is for active variables. 
#' quanti.sup" is for quantitative supplementary variables.
#' @returnType 图形文件
#' @return 
#' 
#' @author dragon 2018.10.7

#' @note 该函数一般用于分析是否数据存在明显的批次效应---观察PCA结果是否直接存在明显的组别差异情况

FunctionPCA <- function(expressionMatrix, groupLabel, plotFile, is.log = TRUE, label = c("all", "none", "ind", "ind.sup", "quali", "var", "quanti.sup")){
        library("FactoMineR")
        library("factoextra")
        if(!is.log){
                expressionMatrix <- log2(expressionMatrix + 1)
        }
        #PCA分析要求行为样本，列为变量
        expressionMatrix <- t(expressionMatrix)
        
        res.pca <- PCA(expressionMatrix, scale.unit = TRUE, graph = F)
        plot1 <- fviz_eig(res.pca)
        # Graph of individuals. Individuals with a similar profile are grouped together.
        plot2 <- fviz_pca_ind(res.pca,
                     col.ind = "cos2", # Color by the quality of representation
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                     repel = TRUE,     # Avoid text overlapping
                     label = label
        )
        if(is.factor(groupLabel)){
                groups <- groupLabel
        }else{
                groups <- as.factor(groupLabel)  
        }

        plot3 <- fviz_pca_ind(res.pca,
                     col.ind = groups, # color by groups
                     addEllipses = TRUE, # Concentration ellipses
                     ellipse.type = "confidence",
                     legend.title = "Groups",
                     repel = TRUE,
                     label = label
        )
        pca<-prcomp(expressionMatrix)
        pca_1 <- predict(pca)
        pca_1 <-as.data.frame(pca_1)
        pca_1$group <-groups
        pca_data_1_prcomp_2<-pca_1[,1:2]# 选取的列数为主成分的个数
        plot4 <- ggplot(pca_data_1_prcomp_2,aes(x=PC1,y=PC2,colour=pca_1$group))+geom_point()
        
        pdf(file = plotFile)
        plot(plot1)
        plot(plot2)
        plot(plot3)
        plot(plot4)
        dev.off()
        
        return(res.pca)
}
