##   TO DO: 芯片表达谱探针转基因表达谱
##
##   Author: jeason zhao
#####################################################

#' @param expr 表达谱，行是探针，列是样本
#' @param ids id转换关系，至少包含两列：ID[探针名]，Gene[基因名]
#' @param type 按中值或者均值合并探针，参数：median，mean

probe2gene <- function(expr, ids, type = c("median", "mean")) {
  library(tidyverse)
  expr %<>% as.data.frame() %>% rownames_to_column(var = "ID_REF")
  meta_expr <- expr %>% 
    dplyr::filter(is.element(ID_REF, ids$ID)) %>% 
    mutate(ID_REF = ids[match(ID_REF, ids$ID),"Gene"]) %>% 
    group_by(ID_REF) %>% nest() %>% ungroup()
  
  type = match.arg(type)
  if (type == "median") {
    meta_expr %<>% mutate(expr = purrr::map(data, ~ apply(.x, 2, median)))
  } else if (type == "mean") {
    meta_expr %<>% mutate(expr = purrr::map(data, ~ apply(.x, 2, mean)))
  } else {
    stop("parameter `type` must be 'median' or 'mean' !")
  }
  res_expr <- do.call(rbind, meta_expr$expr)
  rownames(res_expr) <- meta_expr %>% pull(ID_REF)
  return(res_expr)
}

