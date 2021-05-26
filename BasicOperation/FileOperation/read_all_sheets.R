#########################################
# TODO: 读取一个excel的所有sheets
#
# Author: Jeason Zhao
#########################################

#' @param path Path to the xls/xlsx file.
read_all_sheet <- function(path) {
  if (!requireNamespace("readxl")) {
        stop("Your need install the `readxl` package")
    }
  sheetnames <- readxl::excel_sheets(path)
  res <- lapply(sheetnames, read_excel, path = path)
  names(res) <- sheetnames
  return(res)
}