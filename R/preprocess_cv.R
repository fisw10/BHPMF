PreprocessCv <- function(X, hierarchy.info, num.folds, tmp.dir, verbose) {
  #Missing Values should be stored as NA
  
  # source("countNumTraits.R")
  # source("generatePermutation.R")
  # source("splitData.R")
  # source("findFold.R")
  # source("buildUpperLevelMat.R")
  
  # library(Matrix)
  
  if (!sum(is.na(X))) {
    stop("No missing value in X! Missing Values should be stored as NA.")
  }
  
  if (num.folds < 2) {
    stop("num.folds should be greater than 2.")
  }
  
  if (!is.matrix(X)) {
    stop("X should be a matrix. You may use as.matrix to convert X to a matrix.")
  }
  
  nrows = nrow(X);
  ncols = ncol(X);
  
  if (nrow(hierarchy.info) != nrows) {
    stop("Number of rows of hierarchy.info and X should be the same!")
  }
  
  num.hierarchy.levels <- ncol(hierarchy.info)
  
  column_names <- colnames(X)
  
  # check if all columns has at least one observation. 
  for (ind in 1: ncols) {
    if (!sum(!is.na(X[, ind]))) {
      stop("there is no available observation for column", ind, "!\n")
    }
  }
  
  # check if the hierarchy.info contains all information for all observations
  for (ind in 1: num.hierarchy.levels) {
    if (sum(!is.na(hierarchy.info[, ind])) < nrows) {
      stop("there are observation with no hierarchy info at level", ind, "!\n")
    }
  }
  
  # print("Count Number of Traits per Observation ...")
  # noTraits: a vector of size observation containing number of available trait per observation
  # indObs:  a vector containing index of observation with at least one available trait
  # traitIds: containing index of available trait per observation
  res <- countNumTraits(X, 1);
  noTraits <- res[[1]];
  indObs <- res[[2]];
  traitIds <- res[[3]];
  fileName  <- paste(tmp.dir, "/freqObsNoTraits_pmf.Rda", sep="");
  save(noTraits, indObs ,traitIds, file=fileName)
  
  if (length(indObs) < nrows) {
    stop("There are observations with all features missing! \n 
         Observations with all features missing should be removed!")
  }
  
  # find non-zero elements
  Idx <- which(!is.na(X), arr.ind=TRUE);
  val <- X[Idx];
  
  # Y: 3-column matrix: sparse representation of trait info
  # each row is a triple (i, j, value)
  # i: observation index (row index)
  # j: feature index (column index)
  # vlaue: the observed value
  Y <- cbind(Idx, val);
  rm(Idx, val);
  
  file_name <-paste(tmp.dir, "/ProcessedData_pmf.Rda", sep="")
  save (X, Y, column_names, file = file_name)
  rm(file_name, res)
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # %%Split Dataset to training, test, validation set
  # %% Similar to CV (10 fold) some data may appear more than once in test set
  # %% (not completly disjoint)
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # call generatePermutation
  ptm <- proc.time()
  # permute_file <-paste(output_dir, "/permutationMat.Rda", sep="");
  file_name <-paste(tmp.dir, "/traits_trees_atleastone_", sep="");
  tree_file <-paste(tmp.dir, "/ProcessedData_pmf.Rda", sep="");
  # split the dataset for cross validation
  if(verbose) {
    cat("Splitting the data ...")
  }
  # find_folds(X, num.folds, file_name, tree_file, output_dir, 0, trait_names, 1, num.hierarchy.levels)
  find_folds(X, num.hierarchy.levels, num.folds, 0, tmp.dir)
  
  if (num.hierarchy.levels > 1) {
    if(verbose) {
      cat("preprocess the hierachy: \n")
    }
    # source("findPhyloInfo.R")
    PreprocessHierarchy(X, hierarchy.info, tmp.dir)
    if(verbose) {
      cat("\n Building upper level matrices...")
    }
    # build upper level matrices for each fold
    buildUpperLevelMat(tmp.dir, 0, column_names, num.hierarchy.levels, num.folds, verbose)
    
    if(verbose) {
      cat("preprocessing time: ", proc.time() - ptm, "\n")
    }
  } else {
    num.nodes.per.level <- rep(0, 1)
    num.nodes.per.level[1] <- nrow(X)
    file_name <-paste(tmp.dir, "/processed_hierarchy_info.Rda", sep="")
    save(num.nodes.per.level, file=file_name)
  }
  
  }