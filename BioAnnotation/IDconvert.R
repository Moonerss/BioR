#####################################################
## TODO: 实现基于bioMart和 org.Hs.eg.db的ID注释
## Author： Erjie Zhao
#####################################################

#################################################################################
#' @description 暂时先实现三类ID之间的转换：Ensmbl ID, Symbol, Entrez ID
#' @param ids 字符向量，需要转换的ID
#' @param from 需要转换的ID的类型，可以选：SYMBOL，ENSEMBL，ENTREZID
#' @param to 转换到的ID类型，可以选：SYMBOL，ENSEMBL，ENTREZID
#' @return 返回一个数据框，包含四列，原始ID(origin_id), 处理后的ID（from参数），转换后的ID（to参数）

id_convert <- function(ids, from, to) {
  ## 检测BiocMananger包是否安装
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    stop("There was no `BiocManager` package")
  }
  ## 检测是否安装需要的包
  pkgs <- c("biomaRt", "org.Hs.eg.db", "tidyverse")
  pkgs_index <- is.element(pkgs, installed.packages()[, "Package"])
  if (!all(pkgs_index)) {
    no_pkgs <- pkgs[!pkgs_index]
    if (length(no_pkgs) == 1) {
      answer <- .getAnswer(msg = paste(no_pkgs, "package is not available, install it? [Y/N]", sep = " "),
                           allowed = c("y", "Y", "n", "N"))      
    } else {
      answer <- .getAnswer(msg = paste(paste(no_pkgs, collapse = "、"), "packages is not available, install it? [Y/N]", sep = " "),
                           allowed = c("y", "Y", "n", "N"))
    }
    
    if (answer == "y") {
      install.packages(no_pkgs, dependencies = TRUE)
    } else {
      cat(paste("No package:", paste(no_pkgs, collapse = " ")))
      stop()
    }
  }
  
  ## Step1: 检测ENSG号是否带有版本信息
  if (from == "Ensmbl") {
    if(all(str_detect(ids, pattern = "\\."))) {
      cat("The Ensmbl ID have version, we remove it!")
      ids <- stringr::str_split(rownames(ids),"\\.",simplify = T)[,1]
      gn <- ids[!duplicated(ids)]
      cat(paste("left", length(gn), "Ensmbl ID", collapse = " "))
    }
  }
  
  ## Step2: ID转换，使用两个工具进行
  sapply(pkgs, require, character.only = TRUE)
  from_sclae <- switch(from, "SYMBOL" = "hgnc_symbol", "ENSEMBL" = "ensembl_gene_id", "ENTREZID" = "entrezgene_id")
  
  ## 先用bioMart进行注释,注释不到的就是""
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  convert_id1 <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id"),
                      filters = from_sclae,
                      values = ids,
                      mart = ensembl)
  colnames(convert_id1) <- c("SYMBOL", "ENSEMBL", "ENTREZID")
  
  if (sum(convert_id1[, to] == "") != 0) {
    left_ids <- convert_id1[convert_id1[, to] == "", from]
  }
  convert_id1 <- convert_id1[convert_id1[, to] != "", ]
  
  ## 使用org.Hs.eg.db对未注释到的重新进行注释
  
  all_id <- c("SYMBOL", "ENSEMBL", "ENTREZID")
  convert_id2 <- select(org.Hs.eg.db, keys = left_ids, columns = all_id, keytype = from)
  
  if (sum(convert_id2[, to] == "") != 0) {
    cat(paste(sum(convert_id2[, to] == ""), "ID no match! Keep the origin ID", sep = " "))
  } else {
    cat(paste("All", from, "have been transformed!", sep = " "))
  }
    
  convert_id2[convert_id2[, to] == "", to] <- convert_id2[convert_id2[, to] == "", from]
  
  converted_id <- rbind(convert_id1, convert_id2) 
  
  alls <- bind_cols(ids, converted_id)
  colnames(alls)[1] <- "origin_id"
  
  res <- dplyr::select(alls, c("origin_id", from, to))
  return(res)
}



##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## useful function
##
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
.getAnswer <- function(msg, allowed) {
  if (interactive()) {
    repeat {
      cat(msg)
      answer <- readLines(n = 1)
      if (answer %in% allowed)
        break
    }
    tolower(answer)
  } else {
    "n"
  }
}










