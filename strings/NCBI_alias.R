#######################################################
## TODO: 从ncbi给的基因注释文件种提取基因相对应的名字
## Author: Erjie Zhao
## Time：2020-7-28
## Update：
#######################################################

ncbi_to_alias <- function(filepath, pattern) {
  ncbi <- data.table::fread(filepath, header = T, sep = "\t",stringsAsFactors = F)
  alias <- ncbi$dbXrefs
  split_alias <- strsplit(alias, split = '[|]')
  position <- purrr::map(split_alias, function(x) {
    pos <- stringr::str_which(x, pattern = pattern)
    if (length(pos) > 1) {
      pos <- pos[1]
    }
    pos})
  id_list <- purrr::map2(split_alias, position, function(x,y) {x[y]})
  none <- unlist(lapply(id_list, identical, character(0)))
  id_list[which(none == T)] <- "?"
  ids <- unlist(id_list)
  ids_list <- purrr::map(ids, stringi::stri_split_fixed, pattern = ":", n = 2)
  ids <- unlist(lapply(ids_list, function(x){unlist(x)[2]}))
  ids[which(none == T)] <- "?"
  return(ids)
}
