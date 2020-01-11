##  TO DO: 容易忘记但又经常用到的R操作
##
##  Author: jeason zhao
###########################################################


# 查看已加载的包
(.packages())

# 移除已加载的包
detach("package:xxx")

# 列出包所在库的路径	
.libPaths()

# 更新包
update.packages()

# 卸载包
remove.packages()

# 读取有空值的数据，read.table需要设置
na.strings = ""

# read.table参数，避免读取文件中#字符误当注释处理
comment.char = ""

# 构建固定长度的list
vector("list", n)

# 查看R包的网页版使用介绍  
vignette(package = "packageNames")

# ()内条件是否全部为TRUE
all()

# debug排查函数错误行
options(error = quote({
  #setwd('~/myUsername/directoryForDump'); # Set working directory where you want the dump to go, since dump.frames() doesn't seem to accept absolute file paths.
  dump.frames("errorDump", to.file=TRUE, include.GlobalEnv=TRUE); # First dump to file; this dump is not accessible by the R session.
  sink(file="error.log"); # Specify sink file to redirect all output.
  dump.frames(); # Dump again to be able to retrieve error message and write to error log; this dump is accessible by the R session since not dumped to file.
  cat(attr(last.dump,"error.message")); # Print error message to file, along with simplified stack trace.
  cat('\nTraceback:');
  cat('\n');
  traceback(2); # Print full traceback of function calls with all parameters. The 2 passed to traceback omits the outermost two function calls.
  sink()}))

# 检查是否是list最好使用如下函数（数据框的is.list也是ture）
.is_list <- function(x){
  inherits(x, "list")
}