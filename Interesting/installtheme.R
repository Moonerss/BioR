# RStudio supports user-created Editor themes, including themes in the
# .tmTheme format. Editor themes can be installed from a local file or url 
# address using the {rstudioapi} package. 
#
# Many .tmTheme files can be found in the online repository:
# http://tmtheme-editor.herokuapp.com/#!/editor/theme/Monokai,
# for use in RStudio.  
#
# This file provides options for (1) installing a single theme from a local 
# file, or (2) insttalling all .tmThemes available 

#------------------------------------------------------------------------------#
# Option 1: Install a single theme.
#------------------------------------------------------------------------------#
#   1 - Go to http://tmtheme-editor.herokuapp.com/#!/editor/theme/Monokai
#   2 - Set the language to R.
#   3 - Choose a theme and download it. 
#   4 - Install rstudioapi (if necessary).
#   5 - Install theme using rstudioapi::addTheme(path_to_file)
#   6 - Accecss theme via RStudio > Preview > Appearance > Editor theme (mac)
#------------------------------------------------------------------------------#

#install.packages("rstudioapi")

path_to_file <- "~/Desktop/BBEdit.tmTheme"
rstudioapi::addTheme(path_to_file)

#------------------------------------------------------------------------------#
# Option 2: Install all available themes (takes a few minutes).
#------------------------------------------------------------------------------#
#   1 - Install rstudioapi, jsonlite and httr (if necessary).
#   2 - Run script below.
#   3 - Accecss themes via RStudio > Preview > Appearance > Editor theme (mac)
#------------------------------------------------------------------------------#

#install.packages("jsonlite")
#install.packages("rstudioapi")
#install.packages("httr")

# The urlFileExist function is taken from Patrick's anwser at:
# https://stackoverflow.com/questions/7684771/how-to-check-if-a-file-exists-from-a-url

urlFileExist <- function(url){
  HTTP_STATUS_OK <- 200
  hd <- httr::HEAD(url)
  status <- hd$all_headers[[1]]$status
  list(exists = status == HTTP_STATUS_OK, status = status)
}

# Themes are accessed from the JSON file underlying the tm-Theme Editor (link
# above) developed by Allen Bargi (https://github.com/aziz). 

gallery   <- "https://raw.githubusercontent.com/aziz/tmTheme-Editor/master/app/front/public/gallery.json"
json_file <- jsonlite::fromJSON(txt=gallery,flatten=TRUE)

# Function to install theme from URL with some basic error handling.

installTheme <- function(url,name) {
    out <- tryCatch(
        {
              if(urlFileExist(url)$exists){
                rstudioapi::addTheme(url, force=TRUE)
              }
          
              message(paste("Theme installed:", name))
        },
        error=function(cond) {
            message(paste("Theme caused an error:", name))
            return(NULL)
        },
        warning=function(cond) {
            message(paste("Theme caused a warning:", name))
            return(NULL)
        },
        finally={}
    )    
    return(out)
}

# Install themes.

for (i in 1:nrow(json_file)){
  installTheme(json_file$url[i],json_file$name[i])
}