#' @TODO 列联表检验
#' @param  testTable 列联表
#' @returnType list
#' @return 检验结果，包括 odds ratio, pvalue, 以及计算方法
#'
#' @author ZGX
getChisqTest <- function(testTable=NULL){
	testTable <- as.matrix(testTable)
	# 计算理论频数(theoretical frequency )
	# 理论频数ni不是表格中实际的数，而是计算得到的理论频数。
	# 例如：n1的理论频数即为：(n1+n2)*(n1+n3)/N
	# n2为：(n1+n2)*(n2+n4)/N，依次类推。
	# 计算时，可能会出现数字太大的情况，如下错误
	# NAs produced by integer overflow
	# 因此需要改变运算顺序
	# 计算理论频数函数
	thFreq <- function(o1, o2, total){
		(max(o1, o2)/total)*min(o1, o2)
	}
	# 计算理论频数
	n1 <- thFreq(sum(testTable[1,]), sum(testTable[,1]), sum(testTable))
	n2 <- thFreq(sum(testTable[1,]), sum(testTable[,2]), sum(testTable))
	n3 <- thFreq(sum(testTable[2,]), sum(testTable[,1]), sum(testTable))
	n4 <- thFreq(sum(testTable[2,]), sum(testTable[,2]), sum(testTable))
	theorFreq <- c(n1, n2, n3, n4)
	testRes <- NULL
	# 依据数据，进行不同类别的列联表检验
	if(sum(testTable)>40){
		if(!any(theorFreq<5)){
			# 特别注意：
			# 1）一般情况下使用卡方检验即可（总样本数N=n1+n2+n3+n4> 40且理论频数n1、n2、n3、n4均大于等于5）。R：chisq.test(rbind(c(n1,n2),c(n3,n4)), 2, 2), correct=F)。
			testRes <- chisq.test(testTable, correct=FALSE)
		}
		if(sum(theorFreq<5)==1){
			# 2）若四格表中其中一个数的理论频数ni<5，但N>=40，则使用连续性修正的卡方检验。R：chisq.test(rbind(c(n1,n2),c(n3,n4)), 2, 2), correct=T)。
			testRes <- chisq.test(testTable, correct=TRUE)
		}
		if(sum(theorFreq<5)>1 || any(theorFreq<1)){
			# 3）若两个或两个以上理论频数小于5，N > 40，则可以使用Fisher精确检验（也可以用修正的卡方检验，选择结果好的那个，一般Fisher精确检验会更好些）。如果四格表中个数差异太过明显，两种方法结果都不太好，则可以考虑其他方法，如随机抽取。R：fisher.test(rbind(c(n1,n2),c(n3,n4)), 2, 2))
			# 4）其中一个理论频数小于1，用Fisher精确检验更准确。R：fisher.test(rbind(c(n1,n2),c(n3,n4)), 2, 2))
			testRes <- fisher.test(testTable)
		}
	}else{
		# 4）N<40时，或者其中一个理论频数小于1，用Fisher精确检验更准确。R：fisher.test(rbind(c(n1,n2),c(n3,n4)), 2, 2))
		testRes <- fisher.test(testTable)
	}
	library(fmsb)
	# 计算odds ratio
	oddsratioRes <- oddsratio(testTable)
	# oddsratioRes_estimate <- oddsratioRes$estimate
	# colnames(oddsratioRes_estimate)[1:dim(testTable)[2]] <- colnames(testTable)
	# rownames(oddsratioRes_estimate)[1:dim(testTable)[1]] <- rownames(testTable)
	
	# 建立返回向量
	resList <- list(oddsratioRes$estimate, testRes$p.value, oddsratioRes$conf.int, oddsratioRes$method, testRes$method)
	names(resList) <- c("oddsratio", "pvalue", "oddsratio_95CI", "oddsratio_method", "test_method")
	return(resList)
}



