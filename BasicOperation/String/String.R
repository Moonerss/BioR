# TODO: Add comment
# 
# Author: jeason
###############################################################################


#######################replace substring###################
#chartr(old, new, x)
#
#y <- gsub("(...)", "\\1_", y);#inserts "_" after each triplet
#
#
#
#
#
#######################trim string##########################
#strtrim(c("abcdef", "abcdef", "abcdef"), c(1,5,10));#Trim character strings to specified display widths. 
#
#
#######################format string#######################
#strwrap(x, width = 60, indent = 5); Wrap Character Strings to Format Paragraphs
#
#
#
#######################vectorizes single sequence##########
#substring(singleSeq, 1:nchar(singleSeq), 1:nchar(singleSeq));
#
#
##删除头尾空格
#library(stringr)
#str_trim()


#将第一个字母大写其余小写
firstUpOtherLow <- function(x)
{
	paste(toupper(substr(x,1,1)), tolower(substring(x,2)),sep="")
}