#' Automate ABSOLUTE calling for multiple samples in parallel way
#' 
#' If calling for a sample failed, the error message will be written to
#' `error.log` under result directory.
#' [ABSOLUTE](https://www.nature.com/articles/nbt.2203) is a famous software
#' developed by Broad Institute, however the `RunAbsolute` function is designed for
#' computing one sample each time and set no default values. **DoAbsolute** help
#' user set default parameters according to [ABSOLUTE documentation](http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/ABSOLUTE), provide an uniform interface to
#' input data easily and run RunAbsolute parallelly.
#' 
#' More detail about how to analyze ABSOLUTE results please see [this link](http://software.broadinstitute.org/cancer/software/genepattern/analyzing-absolute-data).
#' @param Seg a `data.frame` or a file (path) contains columns
#' "Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean".
#' @param MAF MAF, default is `NULL`, can provided as `data.frame` or file path.
#' @param sigma.p Provisional value of excess sample level variance used for mode search. Default: 0
#' @param max.sigma.h Maximum value of excess sample level variance (Eq. 6). Default: 0.2
#' @param min.ploidy Minimum ploidy value to consider. Solutions implying lower ploidy values will be discarded. Default: 0.5
#' @param max.ploidy Maximum ploidy value to consider. Solutions implying greater ploidy values will be discarded. Default: 8
#' @param primary.disease Primary disease of the sample. Default: `NA`
#' @param platform one of "SNP_6.0", "Illumina_WES", "SNP_250K_STY". Default: "SNP_6.0"
#' @param copy_num_type The type of copy number to be handled. Either total or allelic. Total is what this package for. Default: "total"
#' @param max.as.seg.count Maximum number of allelic segments. Samples with a higher segment count will be flagged as 'failed'. Default: 1500
#' @param max.non.clonal Maximum genome fraction that may be modeled as non-clonal (subclonal SCNA). Solutions implying greater values will be discarded. Default: 0.05
#' @param max.neg.genome Maximum genome fraction that may be modeled as non-clonal with copy-ratio below that of clonal homozygous deletion. Solutions implying greater values will be discarded. Default: 0.005
#' @param min.mut.af Minimum mutation allelic fraction. Mutations with lower allelic fractions will be filtered out before analysis. Default: 0.1
#' @param min.no.mut Minor allele frequency file, or NULL if one is not available. This specifies the data for somatic point mutations to be used by ABSOLUTE. Default: 5
#' @param nThread number of cores used for computation. Default: 1L
#' @param temp_path directory path used to store tempory files. Default: Absolute subdirectory under `tempdir()`
#' @param out_path directory path used to store result files. Default: work directory
#' @param recover if `TRUE`, recover previous unfinished work.
#' @param auto_review if `TRUE`, use the `ExtractReviewedResults` review the optimal model, else keep the summary result and review by yourselves.
#' @param keepAllResult if `TRUE`, clean all results, otherwise clean result directory and keep most important results. Default: `TRUE`
#' @param keep_temp if `TRUE`, keep temp dir at the end. Default: `TRUE`
#' @param verbose if `TRUE`, print extra info. Default: `FALSE`

quick_absolute <- function(Seg, MAF = NULL, sigma.p = 0, max.sigma.h = 0.2,
                           min.ploidy = 0.5, max.ploidy = 8,
                           primary.disease = NA,
                           platform = c("SNP_6.0", "Illumina_WES", "SNP_250K_STY"),
                           copy_num_type = c("total", "allelic"),
                           max.as.seg.count = 1500,
                           max.non.clonal = 0.05, max.neg.genome = 0.005,
                          min.mut.af = 0.1, min.no.mut = 5, nThread = 1L,
                           temp_path = file.path(tempdir(), "Absolute"),
                           out_path = getwd(), recover = TRUE, auto_review = TRUE, 
                           keepAllResult = TRUE, keep_temp = TRUE, verbose = FALSE) {
  library(data.table)
  library(foreach)
  library(parallel)
  library(doParallel)
  
  ## create output directory
  dir.create2 <- function(x) {
    if (!dir.exists(x)) {
      message("Directory ", x, " does not exists.")
      if (!dir.create(x, recursive = TRUE)) {
        if (dir.create(x)) {
          stop("Failed creating directory!")
        }
      }
    }
  }
  
  if (verbose) cat("=> Setting results directory as", out_path, "\n")
  dir.create2(out_path)
  if (verbose) cat("=> Done!\n")
  
  ## remove pervious file
  if (verbose) cat("=> Checking old log file...\n")
  if (file.exists(file.path(out_path, "error.log"))) {
    unlink(file.path(out_path, "error.log"))
    cat("=> Removed previous error log file.\n")
  }
  if (verbose) cat("=> Done!\n")
  
  ## check the package
  if (verbose) cat("=> Checking the needed package...\n")
  if (!suppressMessages(requireNamespace("ABSOLUTE"))) {
    stop("Find no package called 'ABSOLUTE', please install it...")
  }
  
  if (!suppressMessages(requireNamespace("foreach"))) {
    warning("Find no package called 'foreach', try to install it...")
    install.packages("furrr", dependencies = TRUE)
  }
  
  if (!suppressMessages(requireNamespace("doParallel"))) {
    warning("Find no package called 'doParallel', try to install it...")
    install.packages("parallel", dependencies = TRUE)
  }
  
  if (!suppressMessages(requireNamespace("data.table"))) {
    warning("Find no package called 'data.table', try to install it...")
    install.packages("data.table", dependencies = TRUE)
  }
  if (verbose) cat("=> Done!\n")
  
  ## load segmentation data
  if (verbose) cat("=> Loading segmentation data...\n")
  if (is.character(Seg)) {
    if (file.exists(Seg)) {
      Seg <- data.table::fread(input = Seg)
    } else {
      stop("file ", Seg, " does not exist")
    }
  } else if (inherits(Seg, "data.frame")) {
    Seg <- data.table::setDT(Seg)
  } else {
    Stop("Unsupport Segmentation Format!")
  }
  if (verbose) cat("=> Done!\n")
  
  ## try to load MAF file
  if (!is.null(MAF)) {
    if (verbose) cat("=> Loading Maf data...\n")
    if (is.character(MAF)) {
      if (file.exists(MAF)) {
        MAF <- data.table::fread(input = MAF)
      } else {
        stop("file ", MAF, " does not exist")
      }
    } else if (inherits(MAF, "data.frame")) {
      MAF <- data.table::setDT(MAF)
    } else {
      Stop("Unsupport MAF Format!")
    }
    if (verbose) cat("=> Done!\n")
  }
  
  ## check the data
  if (verbose) cat("=> Checking data format of segmentation file...\n")
  # check segment file
  seg_cols <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
  if (!all(seg_cols %in% colnames(Seg))) {
    stop("Columns ", paste(seg_cols, collapse = " "), " are should in file.")
  } else {
    Seg <- Seg[, c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")]
  }
  
  Seg$Chromosome <- as.character(Seg$Chromosome)
  Seg$Chromosome <- gsub(pattern = "chr", replacement = "", Seg$Chromosome, ignore.case = TRUE)
  Seg$Chromosome <- gsub(pattern = "X", replacement = "23", Seg$Chromosome, ignore.case = TRUE)
  if (verbose) cat("==> Keeping only chr 1-23 for CNV data...\n")
  autosome <- as.character(seq(1, 23))
  Seg <- Seg[Chromosome %in% autosome, ]
  if (verbose) cat("=> Done!\n")
  
  # check maf file
  if (!is.null(MAF)) {
    if (verbose) cat("=> Checking data format of maf file...\n")
    maf_cols <- c(
      "Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_position", "dbSNP_Val_Status",
      "t_ref_count", "t_alt_count"
    )
    if (all(maf_cols %in% colnames(MAF))) {
      MAF <- MAF[, c(
        "Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_position", "dbSNP_Val_Status",
        "t_ref_count", "t_alt_count"
      )]
    } else {
      stop("Necessary columns for Maf file: \n", paste(maf_cols, collapse = " "), "\n")
    }
    
    MAF$Chromosome <- as.character(MAF$Chromosome)
    MAF$Chromosome <- gsub(pattern = "chr", replacement = "", MAF$Chromosome, ignore.case = TRUE)
    MAF$Chromosome <- gsub(pattern = "X", replacement = "23", MAF$Chromosome, ignore.case = TRUE)
    if (verbose) cat("==> Keeping only chr 1-23 for Maf data...\n")
    MAF <- MAF[Chromosome %in% autosome, ]
    if (verbose) cat("=> Done!\n")
  }
  
  ## keep a directory to save temp file
  if (verbose) cat("=> Creating temp directory ...\n")
  dir.create2(temp_path)
  if (verbose) cat("=> Done!\n")
  
  ## split segment data into different files
  if (verbose) cat("=> Spliting seg data of samples to different files...\n")
  samples <- unique(Seg$Sample)
  SAMPLE <- samples
  seg_filepath <- vector(mode = "character", length = length(samples))
  maf_filepath <- vector(mode = "character", length = length(samples))
  
  for (i in seq_along(samples)) {
    if (verbose) cat("==> Processing sample ", samples[i], "...\n")
    seg <- Seg[Sample == samples[i], ]
    seg_filepath[i] <- file.path(temp_path, paste0(samples[i], ".seg"))
    data.table::fwrite(x = seg, file = seg_filepath[i], sep = "\t")
    
    if (is.null(MAF)) {
      maf_filepath[i] <- NA_character_
    } else {
      maf <- MAF[Tumor_Sample_Barcode == samples[i], ]
      if (nrow(maf) == 0) {
        if (verbose) cat("===> This sample has not Maf data, skipping...\n")
        maf_filepath[i] <- NA_character_
      } else {
        if (verbose) cat("===> Filtering mutations which vaf<", min.mut.af, "...\n")
        maf <- maf[(t_alt_count / (t_ref_count + t_alt_count)) >= min.mut.af, ]
        if (verbose) cat("===> Filtering Maf which count<", min.no.mut, "...\n")
        if (nrow(maf) < min.no.mut) {
          if (verbose) cat("===> This sample has not Maf data (after filtering), skipping...\n")
          maf_filepath[i] <- NA_character_
        } else {
          maf_filepath[i] <- file.path(temp_path, paste0(samples[i], ".maf"))
          data.table::fwrite(x = maf, file = maf_filepath[i], sep = "\t")
        }
      }
    }
  }
  if (verbose) cat("=> Done!\n")
  
  ## match options
  platform <- match.arg(platform)
  copy_num_type <- match.arg(copy_num_type)
  
  cache_dir <- file.path(temp_path, "cache")
  dir.create2(cache_dir)
  
  ## get the unfinished work
  if (recover) {
    if (verbose) cat("=> Checking samples have been called...\n")
    not_called <- c()
    maf_filepath <- c()
    seg_filepath <- c()
    for (i in seq_along(samples)) {
      if (!file.exists(file.path(cache_dir, paste0(samples[i], ".ABSOLUTE.RData")))) {
        if (verbose) cat("==> ", samples[i], "is not called.\n")
        not_called <- c(not_called, samples[i])
        if (!file.exists(file.path(temp_path, paste0(samples[i], ".seg")))) {
          stop("Something wrong with splitted files.")
        }
        seg_filepath <- c(seg_filepath, file.path(temp_path, paste0(samples[i], ".seg")))
        if (!file.exists(file.path(temp_path, paste0(samples[i], ".maf")))) {
          maf_filepath <- c(maf_filepath, NA_character_)
        } else {
          maf_filepath <- c(maf_filepath, file.path(temp_path, paste0(samples[i], ".maf")))
        }
      }
    }

    samples <- not_called
    if (length(not_called) != 0) {
      if (verbose) cat("==> ABSOLUTE calling for above samples will be recovered.\n")
    } else {
      if (verbose) cat("==> ABSOLUTE calling has been done. \n")
    }
    if (verbose) cat("=> Done!\n")    
  }
  
  ## run absolute
  if (length(samples) != 0) {
    if (verbose) cat("=> Running RunAbsolute...(by patient)\n")
    nCores <- parallel::detectCores(logical = FALSE)
    if (nThread > nCores) {
      warning("Number of real physical core is ", nCores, " while user set ", nThread, "\nUse ", nCores, " Cores.")
      nThread <- nCores
    }
    
    if (nThread == 1) {
      maf_fn <- ifelse(is.na(maf_filepath), NULL, maf_filepath)
      seg_fn <- seg_filepath
      for (i in seq_along(samples)) {
        if (verbose) cat("==> Processing sample ", samples[i], "...\n")
        tryCatch(
          {
            suppressWarnings(ABSOLUTE::RunAbsolute(
              seg.dat.fn = seg_fn, maf.fn = maf_fn,
              sample.name = samples[i],
              sigma.p = sigma.p, max.sigma.h = max.sigma.h,
              min.ploidy = min.ploidy, max.ploidy = max.ploidy,
              primary.disease = primary.disease, platform = platform,
              results.dir = cache_dir, max.as.seg.count = max.as.seg.count,
              max.non.clonal = max.non.clonal, max.neg.genome = max.neg.genome,
              copy_num_type = copy_num_type,
              min.mut.af = min.mut.af, verbose = verbose
            ))
          },
          error = function(e) {
            cat("===> Error detected, see log for more details.\n")
            sink(file.path(out_path, "error.log"), append = TRUE)
            cat("Detected error in sample", samples[i], "\n")
            cat("Error message:", e$message, "\n")
            if (grepl("mutations left", e$message)) {
              cat("Try fixing by removing Maf file.\n")
              tryCatch(
                {
                  suppressWarnings(ABSOLUTE::RunAbsolute(
                    seg.dat.fn = seg_fn, maf.fn = NULL,
                    sample.name = samples[i],
                    sigma.p = sigma.p, max.sigma.h = max.sigma.h,
                    min.ploidy = min.ploidy, max.ploidy = max.ploidy,
                    primary.disease = primary.disease, platform = platform,
                    results.dir = cache_dir, max.as.seg.count = max.as.seg.count,
                    max.non.clonal = max.non.clonal, max.neg.genome = max.neg.genome,
                    copy_num_type = copy_num_type,
                    min.mut.af = NULL, verbose = verbose
                  ))
                  cat("Fixing successfully!\n")
                },
                error = function(e) {
                  cat("Fixing failed. Skipping this sample.\n")
                }
              )
            } else {
              cat("Skipping this sample.\n")
            }
            cat("========\n")
            sink()
          }
        )
      }
    } else {
      if (Sys.info()[["sysname"]] == "Windows") {
        cl <- parallel::makeCluster(nThread)
        doParallel::registerDoParallel(cl)
      } else {
        doParallel::registerDoParallel(cores = nThread)
      }
      
      foreach(i = seq_along(samples)) %dopar% {
        if (is.na(maf_filepath[i])) {
          maf_fn <- NULL
        } else {
          maf_fn <- maf_filepath[i]
        }
        seg_fn <- seg_filepath[i]
        if (verbose) cat("==> Processing sample ", samples[i], "...\n")
        tryCatch(
          {
            suppressWarnings(ABSOLUTE::RunAbsolute(
              seg.dat.fn = seg_fn, maf.fn = maf_fn,
              sample.name = samples[i],
              sigma.p = sigma.p, max.sigma.h = max.sigma.h,
              min.ploidy = min.ploidy, max.ploidy = max.ploidy,
              primary.disease = primary.disease, platform = platform,
              results.dir = cache_dir, max.as.seg.count = max.as.seg.count,
              max.non.clonal = max.non.clonal, max.neg.genome = max.neg.genome,
              copy_num_type = copy_num_type,
              min.mut.af = min.mut.af, verbose = verbose
            ))
          },
          error = function(e) {
            cat("----> Error detected, see log for more details.\n")
            sink(file.path(out_path, "error.log"), append = TRUE)
            cat("Detected error in sample", samples[i], "\n")
            cat("Error message:", e$message, "\n")
            if (grepl("mutations left", e$message)) {
              cat("Try fixing by removing Maf file.\n")
              tryCatch(
                {
                  suppressWarnings(ABSOLUTE::RunAbsolute(
                    seg.dat.fn = seg_fn, maf.fn = NULL,
                    sample.name = samples[i],
                    sigma.p = sigma.p, max.sigma.h = max.sigma.h,
                    min.ploidy = min.ploidy, max.ploidy = max.ploidy,
                    primary.disease = primary.disease, platform = platform,
                    results.dir = cache_dir, max.as.seg.count = max.as.seg.count,
                    max.non.clonal = max.non.clonal, max.neg.genome = max.neg.genome,
                    copy_num_type = copy_num_type,
                    min.mut.af = NULL, verbose = verbose
                  ))
                  cat("Fixing successfully!\n")
                },
                error = function(e) {
                  cat("Fixing failed. Skipping this sample.\n")
                }
              )
            } else {
              cat("Skipping this sample.\n")
            }
            cat("========\n")
            sink()
          }
        )
      }
    }
    if (verbose) cat("=> Done!\n")
  }
  
  # Create Review Object in absolute
  absolute_files <- file.path(cache_dir, paste0(SAMPLE, ".ABSOLUTE.RData"))
  if (verbose) {
    cat("=> Files in cache directory:\n")
    print(dir(cache_dir))
  }
  ## check whether file exist
  if (verbose) cat("=> Checking result files...\n")
  file_label <- ifelse(file.exists(absolute_files), TRUE, FALSE)
  real_file <- absolute_files[file_label]
  
  if (length(setdiff(absolute_files, real_file)) > 0) {
    warning("==> The next files do not exist:\n", setdiff(absolute_files, real_file))
  }
  
  if (length(real_file) < 1) {
    stop("No result file to proceed.")
  }
  
  if (verbose) cat("=> Done!\n")
  
  ## remove previous review directory 
  review_dir <- file.path(cache_dir, "review")
  # if (dir.exists(review_dir)) {
  #   if (verbose) cat("=> Removed previous temp review result directory.\n")
  #   unlink(review_dir, recursive = TRUE)
  # }
  # if (verbose) cat("=> Done!\n")
  
  ## Absolute summarize
  if (verbose) cat("=> Running Absolute summarize...\n")
  suppressWarnings(ABSOLUTE::CreateReviewObject(
    obj.name = "quick_absolute",
    absolute.files = real_file,
    indv.results.dir = review_dir,
    copy_num_type = copy_num_type,
    plot.modes = TRUE, verbose = verbose
  ))
  if (verbose) cat("\n=> Absolute summarize done!")

  # whether auto review models
  if (auto_review) {
    if (verbose) cat("=> auto_review is TRUE, prepare auto-reviewing...\n")
    pp_call_fn <- file.path(review_dir, "quick_absolute.PP-calls_tab.txt")
    modes_fn <- file.path(review_dir, "quick_absolute.PP-modes.data.RData")
    suppressWarnings(ABSOLUTE::ExtractReviewedResults(
      reviewed.pp.calls.fn = pp_call_fn,
      analyst.id = "wsx",
      modes.fn = modes_fn,
      out.dir.base = review_dir,
      obj.name = "quick_absolute",
      copy_num_type = copy_num_type,
      verbose = verbose
    ))
    
    if (verbose) cat("=> Done!\n")
  }
  
  # output final results
  reviewed_dir <- file.path(review_dir, "reviewed")
  if (auto_review) {
    if (verbose) cat("=> Outputing final results...\n")
  } else {
    if (verbose) cat("\n=> Outputing final results...\n")
  }
  
  ## whether keep all result
  if (keepAllResult) {
    if (verbose) cat("==> keepAllResult is TRUE, keeping all results...\n")
  }
  
  if (verbose) cat("==> Copying files in temp directory to result directory...\n")
  
  if (auto_review) {
    ## save result
    seg_dir <- file.path(out_path, "seg")
    maf_dir <- file.path(out_path, "maf")
    dir.create2(seg_dir)
    dir.create2(maf_dir)
    turn_file <- file.copy(file.path(reviewed_dir, grep("quick_absolute", dir(reviewed_dir), value = TRUE)), out_path)
    
    
    files_seg <- file.path(
      file.path(reviewed_dir, "SEG_MAF"),
      grep("segtab.txt", dir(file.path(reviewed_dir, "SEG_MAF")), value = TRUE)
    )
    if (verbose) cat("===> Copying ", paste(files_seg, collapse = ", "), "to", seg_dir, "...\n")
    turn_file <- file.copy(from = files_seg, to = seg_dir)
    
    files_maf <- file.path(
      file.path(reviewed_dir, "SEG_MAF"),
      grep("ABS_MAF.txt", dir(file.path(reviewed_dir, "SEG_MAF")), value = TRUE)
    )
    if (verbose) cat("===> Copying ", paste(files_maf, collapse = ", "), "to", maf_dir, "...\n")
    turn_file <- file.copy(from = files_maf, to = maf_dir)
    
    if (keepAllResult) {
      call_dir <- file.path(out_path, "summary")
      samples_before_call_dir <- file.path(out_path, "sample_before_summary")
      samples_called_dir <- file.path(out_path, "sample_final_called")
      dir.create2(call_dir)
      dir.create2(samples_before_call_dir)
      dir.create2(samples_called_dir)
      
      files_call <- file.path(review_dir, grep("quick_absolute", dir(review_dir), value = TRUE))
      turn_file <- file.copy(from = files_call, to = call_dir)
      
      files_before <- file.path(cache_dir, grep("ABSOLUTE", dir(cache_dir), value = TRUE))
      turn_file <- file.copy(from = files_before, to = samples_before_call_dir)
      
      files_samplesCalled <- file.path(
        file.path(reviewed_dir, "samples"),
        dir(file.path(reviewed_dir, "samples"))
      )
      turn_file <- file.copy(from = files_samplesCalled, samples_called_dir)
    } else {
      samples_called_dir <- file.path(out_path, "sample_final_called")
      dir.create2(samples_called_dir)
      
      files_samplesCalled <- file.path(
        file.path(reviewed_dir, "samples"),
        dir(file.path(reviewed_dir, "samples"))
      )
      turn_file <- file.copy(from = files_samplesCalled, samples_called_dir)
    }
  } else {
    if (keepAllResult) {
      call_dir <- file.path(out_path, "summary")
      samples_before_call_dir <- file.path(out_path, "sample_before_summary")
      dir.create2(call_dir)
      dir.create2(samples_before_call_dir)
      
      files_call <- file.path(review_dir, grep("quick_absolute", dir(review_dir), value = TRUE))
      turn_file <- file.copy(from = files_call, to = call_dir)
      
      files_before <- file.path(cache_dir, grep("ABSOLUTE", dir(cache_dir), value = TRUE))
      turn_file <- file.copy(from = files_before, to = samples_before_call_dir)
    } else {
      call_dir <- file.path(out_path, "summary")
      dir.create2(call_dir)
      
      files_call <- file.path(review_dir, grep("quick_absolute", dir(review_dir), value = TRUE))
      turn_file <- file.copy(from = files_call, to = call_dir)
    }
  }
  
  if (verbose) cat("=> Done!\n")
  
  # remove temp directory
  if (!keep_temp) {
    unlink(temp_path, recursive = TRUE, force = TRUE)
  }
  cat("================> ALL is Done<==================")
}
