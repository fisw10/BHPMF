CalculateCvRmse <- function(X, hierarchy.info, num.folds=10, prediction.level,
                            used.num.hierarchy.levels, num.samples=1000, burn=100, gaps=2,
                            num.latent.feats=15, tuning=FALSE, num.folds.tuning, tmp.dir, verbose=FALSE) {
  # This function calculate the average RMSE in cross validation.
  #
  # Args:
  #   X: A matrix containing the missing values, rows are the observations,
  #      and columns are the features. Missing values are filled by NA.
  #   hierarchy.info: A matrix containing the hierarchical info
  #      in the format of taxa-table.
  #   num.folds: Number of cross validation (CV) folds. The default value is 10.
  #   prediction.level: The level at which gaps are filled.
  #      e.g., =4, filling gaps at leaves with a hierarchy informaiton of 3 levels.
  #		 The default value is at the observation level.
  #   used.num.hierarchy.levels: Number of hierarchy levels that is used for gap filling.
  #      e.g., =2 means using only the first and second level of hierarchy.
  #      It should between 1 and total number of hierarchy.
  #		 The default value is total number of hierarchy level.
  #   num.samples: Total number of generated samples at each fold using gibbs sampling.
  #      It is not the effective number of samples.
  #		 The default value is 1000.
  #   burn: Number of initial sampled parameters discarded. The default value is 100.
  #   gaps: Gap between sampled parameters kept. The default value is 2.
  #	  num.latent.feats: Size of latent vectors in BHPMF. Set it if tuning is False.
  #	  	 The default value is 10.
  #	  tuning: If set true, first tune BHPMF to choose the best value of num.latent.feats.
  #	  	 The default value is False.
  #   num.folds.tuning: Number of cross validation folds for tuning. The default value is
  #     10 if tuning set to TRUE.
  #	  tmp.dir: A temporary directory used to save preprocessing files.
  #	  	 If not provided, a tmp directory in /R/tmp will be created.
  #		 If provided, each time calling a function from this package, will use the saved
  #		 preprocessing files saved in this directory, it helps to avoid running preprocessing
  #		 functions for the same input.
  #		 WARNING: This directory should be empty for the first time call on a new dataset.
  #    verbose: if TRUE, progress of the sampler is printed to the screen. 
  #       Otherwise, nothing is printed to the screen. Default is FALSE.
  #
  # Returns:
  #   A list of the average RMSE along with standard deviation in cross validation.
  #   If tuning is true, it also returns the tuned number of latent factors paramter
  #   for BHMPF along with minimum RMSE.
  
  # source("file_exists.R")
  # source("preprocess_cv.R")
  num.nodes.per.level <- NULL
  num.hierarchy.levels <- ncol(hierarchy.info)
  num.cols <- ncol(X)
  
  # Check missing arguments
  if (missing(X)) {
    stop("Missing X!")
  }
  if (missing(hierarchy.info)) {
    stop("Missing hierarchy.info!")
  }
  if (missing(used.num.hierarchy.levels)) {
    used.num.hierarchy.levels <- num.hierarchy.levels - 1
  }
  if (missing(prediction.level)) {
    prediction.level <- num.hierarchy.levels
  }
  if (burn > num.samples) {
    stop("burn should be smaller than number of samples.")
  }
  if (gaps > num.samples) {
    stop("gaps should be smaller than number of samples.")
  }
  preprocess.flag <- FALSE
  if (missing(tmp.dir)) {
    tmp.dir <- paste(tempdir(), "/BHPMFAuthorsTmp", sep = "")
    if (file.exists(tmp.dir)) {
      unlink(tmp.dir, recursive = TRUE, force = TRUE)
    }
    dir.create(tmp.dir)
    preprocess.flag <- TRUE
  } else if (!file.exists(tmp.dir)) {
    stop("tmp.dir: ", tmp.dir, "  does not exist")
  }
  
  if (!preprocess.flag) {   # tmp directory is provided by user
    if (!CheckPreprocessFilesExist(tmp.dir, num.folds, num.hierarchy.levels, prediction.level)) { # return TRUE if files exists
      preprocess.flag <- TRUE
    }
  }
  
  if (preprocess.flag) {
    PreprocessCv(X, hierarchy.info, num.folds, tmp.dir, verbose)
    tmp.tune.dir = paste(tmp.dir, "/fold1/Tunning", sep="")
    dir.create(tmp.tune.dir)
  }
  
  load(paste(tmp.dir, "/processed_hierarchy_info.Rda", sep=""))
  
  #Tune the num.latent.feats parameter
  #By choosing the best parameter that minimize the
  #CV RMSE over training data in the first fold
  #Assumption: the best parameter is the same for other folds
  if (tuning) {
    if (verbose) {
      cat("tuning: ")
    }
    # find the X matrix for fold 1
    file.name <- paste(tmp.dir, '/fold1/Ytrain', prediction.level, ".txt", sep = "")
    Y.tune <- as.matrix(read.table(file.name, sep="\t", header=F))
    # X.tune <- sparseMatrix (Y.tune[,1], Y.tune[,2], x=Y.tune[,3])
    X.tune <- matrix(data=NA, nrow=nrow(X), ncol=ncol(X))
    X.tune[cbind(Y.tune[, 1], Y.tune[, 2])] <- Y.tune[, 3]
    rm(Y.tune)
    
    tmp.tune.dir = paste(tmp.dir, "/fold1/Tunning", sep="")
    
    out.tuning <- TuneBhpmf(X.tune, hierarchy.info, prediction.level, used.num.hierarchy.levels,
                     num.folds.tuning, num.samples, burn, gaps, tmp.tune.dir, verbose)
    rm(X.tune)
    num.latent.feats <- out.tuning$BestNumLatentFeats
  }
  
  rmse_vec = rep(0, num.folds);
  save.file.flag <- 0
  out.whole.flag <- 0
  for (fold in 1 : num.folds) {
    tmp.env <- new.env()
    args <- list(
      "NumSamples" = as.integer(num.samples),
      "InputDir" = tmp.dir,
      "DatasetId" = as.integer(fold),
      "Gaps" = as.integer(gaps),
      "Burn" = as.integer(burn),
      "SaveFileFlag" = as.integer(save.file.flag),
      "OutWholeFlag" = as.integer(out.whole.flag),
      "NumTraits" = as.integer(num.cols),
      "NumFeats" = as.integer(num.latent.feats),
      "NumHierarchyLevel" = as.integer(num.hierarchy.levels-1),
      "PredictLevel" = as.integer(prediction.level),
      "UsedNumHierarchyLevel" = as.integer(used.num.hierarchy.levels),
      "Verbose" = as.integer(verbose),
      "Env" = tmp.env
    )
    
    out <- .Call("DemoHPMF", args, "", "", num.nodes.per.level);
    
    if (!save.file.flag) {
      rmse_vec[fold] <- out$RMSE
      next;
    }
  }
  
  avg_rmse = mean(rmse_vec)
  std_rmse = sd(rmse_vec)
  
  if (tuning) {
    result <- list("avg.rmse" = as.numeric(avg_rmse),
                 "std.rmse" = as.numeric(std_rmse),
                 "min.rmse" = out.tuning$min.rmse,
                 "best.number.latent.features" = out.tuning$mbest.number.latent.features
              )
  } else {
    result <- list("avg.rmse" = as.numeric(avg_rmse),
                   "std.rmse" = as.numeric(std_rmse)
              )
  }
  return(result);
}