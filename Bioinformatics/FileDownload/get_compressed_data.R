#' Download and decompress .zip, .tar, .gz, and .tar.gz files
#'
#' This function to download and decompress a file to a target directory
#' @param url The URL to the .zip file resource
#' @param target_dir The directory where the .zip file should be extracted
#' @param format the format of compress files, support: tar, gz, zip; if don't know the format, choose NULL.
#' @param force Should the data be downloaded even if it exists on disk? Default is FALSE.
#' @keywords tar
#' @export
#' @examples
#' \dontrun{
#' # Download/ extract the `switchboard.zip` file
#' get_compressed_data(url = "https://raw.githubusercontent.com/nltk/nltk_data/gh-pages/packages/corpora/switchboard.zip",
#' target_dir = "./data/scs/")
#'
#' # Download/ extract the `swb1_dialogact_annot.tar.gz` file
#' get_compressed_data(url = "https://catalog.ldc.upenn.edu/docs/LDC97S62/swb1_dialogact_annot.tar.gz",
#' target_dir = "./data/scd/")
#' }
#' get_compressed_data()

get_compressed_data <- function(url, target_dir, format = c("zip", "gz", "tar"), force = FALSE) {
  # Get the extension of the target file
  ext <- tools::file_ext(url)
  if (ext == "") {
    if (!is.null(format)) {
      ext <- match.arg(format)
    }
  }
  # Check to see if the target file is a compressed file
  if (!ext %in% c("zip", "gz", "tar")) stop("Target file given is not supported")
  # Check to see if the data already exists
  if (!dir.exists(target_dir) | force == TRUE) { # if data does not exist, download/ decompress
    cat("Creating target data directory \n") # print status message
    dir.create(path = target_dir, recursive = TRUE, showWarnings = FALSE) # create target data directory
    cat("Downloading data... \n") # print status message
    temp <- tempfile() # create a temporary space for the file to be written to
    download.file(url = url, destfile = temp) # download the data to the temp file
    # Decompress the temp file in the target directory
    if (ext == "zip") {
      unzip(zipfile = temp, exdir = target_dir, junkpaths = TRUE) # zip files
    } else {
      untar(tarfile = temp, exdir = target_dir) # tar, gz files
    }
    cat("Data downloaded! \n") # print status message
  } else { # if data exists, don't download it again
    cat("Data already exists \n") # print status message
  }
}
